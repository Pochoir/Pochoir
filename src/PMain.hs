{-
 ----------------------------------------------------------------------------------
 -  Copyright (C) 2010-2011  Massachusetts Institute of Technology
 -  Copyright (C) 2010-2011  Yuan Tang <yuantang@csail.mit.edu>
 - 		                     Charles E. Leiserson <cel@mit.edu>
 - 	 
 -   This program is free software: you can redistribute it and/or modify
 -   it under the terms of the GNU General Public License as published by
 -   the Free Software Foundation, either version 3 of the License, or
 -   (at your option) any later version.
 -
 -   This program is distributed in the hope that it will be useful,
 -   but WITHOUT ANY WARRANTY; without even the implied warranty of
 -   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 -   GNU General Public License for more details.
 -
 -   You should have received a copy of the GNU General Public License
 -   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 -
 -   Suggestsions:                  yuantang@csail.mit.edu
 -   Bugs:                          yuantang@csail.mit.edu
 -
 --------------------------------------------------------------------------------
 -}

module Main where

import Prelude hiding (catch)
import System.Cmd
import System.IO
import System.FilePath
import System.Environment
import System.Exit
import qualified Control.Exception as Exception
import Data.List
import System.Directory 
import System.Cmd (rawSystem)
import Data.Char (isSpace)
import qualified Data.Map as Map
import Text.ParserCombinators.Parsec (runParser)

import PUtils
import PData
import PMainParser

main :: IO ()
main = do args <- getArgs
          cc <- getEnv "CC_CILK"

          whilst (null args) $ do
             printUsage
             exitFailure
          let (inFiles, inDirs, mode, showFile, userArgs) 
                = parseArgs ([], [], PDefault, True, []) args
          whilst (mode == PHelp) $ do
             printOptions
             exitFailure
          putStrLn ("Before ppopp")
          whilst (mode /= PNoPP) $ do
             ppopp (mode, showFile, userArgs) (zip inFiles inDirs)
          -- pass everything to g++ after preprocessing and Pochoir optimization
          putStrLn ("After ppopp")
          let icc = if cc == "icc" then cc_icc else cc_gcc
          let ccFlags = 
                if cc == "icc"  
                   then if mode == PDebug
                           then debugCCFlags_icc 
                           else ccFlags_icc
                   else if mode == PDebug
                           then debugCCFlags_gcc
                           else ccFlags_gcc
          let iccArgs = ccFlags ++ userArgs
          putStrLn (icc ++ " " ++ intercalate " " iccArgs)
          rawSystem icc iccArgs
          whilst (showFile == False) $ do
             let outFiles = map (rename "_pochoir") inFiles 
             let kernelFiles = map (rename "_kernel_info") inFiles
             removeFile $ intercalate " " outFiles
             removeFile $ intercalate " " kernelFiles

whilst :: Bool -> IO () -> IO ()
whilst True action = action
whilst False action = return () 

ppopp :: (PMode, Bool, [String]) -> [(String, String)] -> IO ()
ppopp (_, _, _) [] = return ()
ppopp (mode, showFile, userArgs) ((inFile, inDir):files) = 
    do putStrLn ("pochoir called with mode =" ++ show mode)
       -- System use environment variables
       cc <- getEnv "CC_CILK"

       whilst (cc == "EnvError") $ do
          putStrLn ("Pochoir environment variable not set correctly!")
          exitFailure
       -- let envPath = ["-I" ++ pochoirLibPath] ++ ["-I" ++ cilkHeaderPath]
       let iccPPFile = inDir ++ getPPFile inFile
       let ppFile = ["-o", iccPPFile]
       let outFile = inDir ++ getPochoirFile inFile
       let kernelFile = rename "_kernel_info" inFile
       let ppFlags = if cc == "icc" 
                        then if mode == PDebug
                                then debugPPFlags_icc ++ [inFile] ++ ppFile
                                else ppFlags_icc ++ [inFile] ++ ppFile
                        else if mode == PDebug 
                                then debugPPFlags_gcc ++ [inFile] ++ ppFile
                                else ppFlags_gcc ++ [inFile] ++ ppFile
       -- a pass of c++ preprocessing
       let cmd = cc ++ " " ++ intercalate " " ppFlags
       putStrLn cmd
       system cmd
       -- a pass of pochoir compilation
       whilst (mode /= PDebug) $ do
           inh <- openFile iccPPFile ReadMode
           outh <- openFile outFile WriteMode
           kernelh <- openFile kernelFile WriteMode
           putStrLn ("pochoir " ++ show mode ++ " " ++ iccPPFile)
           pProcess mode (inFile, inDir) inh outh kernelh
           hClose inh
           hClose outh
           hClose kernelh
       ppopp (mode, showFile, userArgs) files

getPPFile :: String -> String
getPPFile fname = name ++ ".i"
    where (name, suffix) = break ('.' ==) fname

getPochoirFile :: String -> String
getPochoirFile fname = name ++ "_pochoir.cpp"
    where (name, suffix) = break ('.' ==) fname

getBinFile :: String -> String
getBinFile fname = name 
    where (name, suffix) = break ('.' ==) fname

pInitState = ParserState { pMode = PCaching, pInFile = "", pInDir = "", pMacro = Map.empty, pArray = Map.empty, pStencil = Map.empty, pShape = Map.empty, pRange = Map.empty, pKernel = Map.empty, pKernelFunc = Map.empty, pGuard = Map.empty, pGuardFunc = Map.empty, pTile = Map.empty, pTileOrder = 0, pColorNum = 0 }


-- When the environment variable "CC_CILK" equals to "icc", we will use cc_icc; otherwise cc_gcc
cc_icc = "icc"
cc_gcc = "g++_4_7"


-- preprocessing flags for g++
ppFlags_gcc = ["-E", "-DNCHECK_SHAPE", "-DNDEBUG", "-std=c++11", "-Wall", "-fcilkplus", "-I${POCHOIR_LIB_PATH}", "-I${CILK_HEADER_PATH}"]
-- preprocessing flags for icc
ppFlags_icc = ["-P", "-C", "-DNCHECK_SHAPE", "-DNDEBUG", "-std=c++11", "-Wall", "-I/usr/include/x86_64-linux-gnu/", "-I${POCHOIR_LIB_PATH}", "-I${CILK_HEADER_PATH}"]

-- compile flags for gcc
ccFlags_gcc = ["-O3", "-Wall", "-unroll-aggressive", "-funroll-loops", "-std=c++11", "-fcilkplus", "-lcilkrts"]

-- compile flags for icc
ccFlags_icc = ["-DNCHECK_SHAPE", "-DNDEBUG", "-O3", "-Wall", "-unroll-aggressive", "-funroll-loops", "-xHOST", "-fno-alias", "-fno-fnalias", "-std=c++11"]

-- debug flags for icc
debugPPFlags_icc = ["-P", "-C", "-DCHECK_SHAPE", "-DDEBUG", "-std=c++11", "-g", "-Wall", "-I/usr/include/x86_64-linux-gnu/", "-I${POCHOIR_LIB_PATH}", "-I${CILK_HEADER_PATH}", "-include", "cilk_stub.h"]

-- debug flags for g++
debugPPFlags_gcc = ["-E", "-DCHECK_SHAPE", "-DDEBUG", "-std=c++11", "-g", "-Wall", "-I${POCHOIR_LIB_PATH}", "-I${CILK_HEADER_PATH}", "-include", "cilk_stub.h"]

-- debug compile flags for gcc
debugCCFlags_gcc = ["-O0", "-g", "-Wall", "-std=c++11", "-fcilkplus", "-lcilkrts"]

-- debug compile flags for icc
debugCCFlags_icc = ["-O0", "-g", "-Wall", "-std=c++11"]


parseArgs :: ([String], [String], PMode, Bool, [String]) -> [String] -> ([String], [String], PMode, Bool, [String])
parseArgs (inFiles, inDirs, mode, showFile, userArgs) aL 
    | elem "--help" aL =
        let l_mode = PHelp
            aL' = delete "--help" aL
        in  (inFiles, inDirs, l_mode, showFile, aL')
    | elem "-h" aL =
        let l_mode = PHelp
            aL' = delete "-h" aL
        in  (inFiles, inDirs, l_mode, showFile, aL')
    | elem "-auto-optimize" aL =
        let l_mode = PDefault
            aL' = delete "-auto-optimize" aL
        in  (inFiles, inDirs, l_mode, showFile, aL')
------------------------------------------------------------------------------
-- so far, the split-caching mode doesn't work for multiple-kernel case!!!! --
------------------------------------------------------------------------------
    | elem "-split-caching" aL =
        let l_mode = PCaching
            aL' = delete "-split-caching" aL
        in  parseArgs (inFiles, inDirs, l_mode, showFile, aL') aL'
    | elem "-cond-macro" aL =
        let l_mode = PCondMacro
            aL' = delete "-cond-macro" aL
        in  parseArgs (inFiles, inDirs, l_mode, showFile, aL') aL'
    | elem "-cond-c-pointer" aL =
        let l_mode = PCondCPointer
            aL' = delete "-cond-c-pointer" aL
        in  parseArgs (inFiles, inDirs, l_mode, showFile, aL') aL'
    | elem "-cond-pointer" aL =
        let l_mode = PCondPointer
            aL' = delete "-cond-pointer" aL
        in  parseArgs (inFiles, inDirs, l_mode, showFile, aL') aL'
    | elem "-cond-opt-pointer" aL =
        let l_mode = PCondOptPointer
            aL' = delete "-cond-opt-pointer" aL
        in  parseArgs (inFiles, inDirs, l_mode, showFile, aL') aL'
    | elem "-unroll-macro" aL =
        let l_mode = PUnrollMacro
            aL' = delete "-unroll-macro" aL
        in  parseArgs (inFiles, inDirs, l_mode, showFile, aL') aL'
    | elem "-unroll-c-pointer" aL =
        let l_mode = PUnrollCPointer
            aL' = delete "-unroll-c-pointer" aL
        in  parseArgs (inFiles, inDirs, l_mode, showFile, aL') aL'
    | elem "-unroll-pointer" aL =
        let l_mode = PUnrollPointer
            aL' = delete "-unroll-pointer" aL
        in  parseArgs (inFiles, inDirs, l_mode, showFile, aL') aL'
    | elem "-unroll-opt-pointer" aL =
        let l_mode = PUnrollOptPointer
            aL' = delete "-unroll-opt-pointer" aL
        in  parseArgs (inFiles, inDirs, l_mode, showFile, aL') aL'
    | elem "-unroll-multi-kernel" aL =
        let l_mode = PMUnroll
            aL' = delete "-unroll-multi-kernel" aL
        in  parseArgs (inFiles, inDirs, l_mode, showFile, aL') aL'
    | elem "-showFile" aL =
        let l_showFile = True
            aL' = delete "-showFile" aL
        in  parseArgs (inFiles, inDirs, mode, l_showFile, aL') aL'
    | elem "-debug" aL =
        let l_mode = PDebug
            aL' = delete "-debug" aL
        in  parseArgs (inFiles, inDirs, l_mode, showFile, aL') aL'
    | null aL == False =
        let (l_files, l_dirs, l_mode, aL') = findCPP aL ([], [], mode, aL)
        in  (l_files, l_dirs, l_mode, showFile, aL')
    | otherwise = 
        let l_mode = PNoPP
        in  (inFiles, inDirs, l_mode, showFile, aL)

findCPP :: [String] -> ([String], [String], PMode, [String]) -> ([String], [String], PMode, [String])
findCPP [] (l_files, l_dirs, l_mode, l_al) = 
    let l_mode' = 
            if null l_files == True || null l_dirs == True then PNoPP else l_mode
    in  (l_files, l_dirs, l_mode', l_al)
findCPP (a:as) (l_files, l_dirs, l_mode, l_al)
    | isSuffixOf ".cpp" a || isSuffixOf ".cxx" a = 
        let l_file = drop (1 + (pLast $ findIndices (== '/') a)) a
            l_dir  = take (1 + (pLast $ findIndices (== '/') a)) a
            l_files' = l_files ++ [l_file]
            l_dirs'  = l_dirs ++ [l_dir]
            pLast [] = -1
            pLast aL@(a:as) = last aL
            l_pochoir_file = rename "_pochoir" l_file 
            (prefix, suffix) = break (a == ) l_al
            l_al' = prefix ++ [l_pochoir_file] ++ tail suffix
        in  findCPP as (l_files', l_dirs', l_mode, l_al')
    | otherwise = findCPP as (l_files, l_dirs, l_mode, l_al)

printUsage :: IO ()
printUsage =
    do putStrLn ("Usage: pochoir [OPTION] [filename]")
       putStrLn ("Try `pochoir --help' for more options.")

printOptions :: IO ()
printOptions = 
    do putStrLn ("Usage: pochoir [OPTION] [filename]")
       putStrLn ("Run the Pochoir stencil compiler on [filename].")
       putStrLn ("-auto-optimize : " ++ breakline ++ "Let the Pochoir compiler automatically choose the best optimizing level for you! (default)")
       putStrLn ("-split-macro-shadow $filename : " ++ breakline ++ 
               "using macro tricks to split the interior and boundary regions")
       putStrLn ("-split-pointer $filename : " ++ breakline ++ 
               "Default Mode : split the interior and boundary region, and using C-style pointer to optimize the base case")
       putStrLn ("\n--Following modes are for experimental tiled specification of irregular stencil computation--\n")
       putStrLn ("-all-cond-tile-macro $filename : " ++ breakline ++
               "put all conditional check in inner-most loop, and use -split-macro-shadow mode for optimizing the base case")
       putStrLn ("-all-cond-tile-pointer $filename : " ++ breakline ++
               "put all conditional check in inner-most loop, and use -split-pointer mode for optimizing the base case")
       putStrLn ("-unroll-t-tile-macro $filename : " ++ breakline ++
               "unroll the tiled kernels along time dimension, for conditional check on spatial dimension, we leave them in the inner-most kernel, and use -split-macro-shadow mode for optimizing the base case")
       putStrLn ("-unroll-t-tile-pointer $filename : " ++ breakline ++
               "unroll the tiled kernels along time dimension, for conditional check on spatial dimension, we leave them in the inner-most kernel, and use -split-macro-pointer mode for optimizing the base case")

pProcess :: PMode -> (String, String) -> Handle -> Handle -> Handle -> IO ()
pProcess mode (inFile, inDir) inh outh kernelh = 
    do ls <- hGetContents inh
       let pRevInitState = pInitState { pMode = mode, pInFile = inFile, pInDir = inDir }
       case runParser pParser pRevInitState "" $ stripWhite ls of
           Left err -> print err
           -- outh is *_pochoir.cpp; kernelh is *_kernel_info.cpp
           Right (outContent, kernContent) -> 
                do hPutStrLn outh outContent
                   hPutStrLn kernelh kernContent


