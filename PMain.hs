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

import System
import IO hiding (try) -- "try" is also defined in Parsec
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
          whilst (null args) $ do
             printUsage
             exitFailure
          let (inFiles, inDirs, mode, debug, showFile, userArgs) 
                = parseArgs ([], [], PDefault, False, True, []) args
          whilst (mode == PHelp) $ do
             printOptions
             exitFailure
          whilst (mode /= PNoPP) $ do
             ppopp (mode, debug, showFile, userArgs) (zip inFiles inDirs)
          -- pass everything to icc after preprocessing and Pochoir optimization
          let iccArgs = userArgs
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

ppopp :: (PMode, Bool, Bool, [String]) -> [(String, String)] -> IO ()
ppopp (_, _, _, _) [] = return ()
ppopp (mode, debug, showFile, userArgs) ((inFile, inDir):files) = 
    do putStrLn ("pochoir called with mode =" ++ show mode)
       pochoirLibPath <- catch (getEnv "POCHOIR_LIB_PATH")(\e -> return "EnvError")
       whilst (pochoirLibPath == "EnvError") $ do
          putStrLn ("Pochoir environment variable not set:")
          putStrLn ("POCHOIR_LIB_PATH")
          exitFailure
       -- let envPath = ["-I" ++ pochoirLibPath] ++ ["-I" ++ cilkHeaderPath]
       let envPath = ["-I" ++ pochoirLibPath]
       let iccPPFile = inDir ++ getPPFile inFile
       let midFile = getMidFile inFile
       let outFile = rename "_pochoir" inFile
       let kernelFile = rename "_kernel_info" inFile
       let iccPPArgs = if debug == False
             -- then iccPPFlags ++ envPath ++ [inFile] ++ [">", midFile]
             -- else iccDebugPPFlags ++ envPath ++ [inFile] ++ [">", midFile]
             then iccPPFlags ++ envPath ++ [inFile]
             else iccDebugPPFlags ++ envPath ++ [inFile]
       -- a pass of icc preprocessing
       -- rawSystem icc iccPPArgs
       let cmd = icc ++ " " ++ intercalate " " iccPPArgs
       putStrLn cmd
       system cmd
       -- a pass of pochoir compilation
       whilst (mode /= PDebug) $ do
           -- inh <- openFile iccPPFile ReadMode
           inh <- openFile midFile ReadMode
           outh <- openFile outFile WriteMode
           kernelh <- openFile kernelFile WriteMode
           putStrLn ("pochoir " ++ show mode ++ " " ++ iccPPFile)
           pProcess mode (inFile, inDir) inh outh kernelh
           hClose inh
           hClose outh
           hClose kernelh
       whilst (mode == PDebug) $ do
           putStrLn ("mv " ++ midFile ++ " " ++ outFile)
           renameFile midFile outFile
       ppopp (mode, debug, showFile, userArgs) files

getMidFile :: String -> String
getMidFile a  
    | isSuffixOf ".cpp" a || isSuffixOf ".cxx" a = take (length a - 4) a ++ ".i"
    | otherwise = a

getPPFile :: String -> String
getPPFile fname = name ++ ".i"
    where (name, suffix) = break ('.' ==) fname

pInitState = ParserState { pMode = PCaching, pInFile = "", pInDir = "", pMacro = Map.empty, pArray = Map.empty, pStencil = Map.empty, pShape = Map.empty, pRange = Map.empty, pKernel = Map.empty, pKernelFunc = Map.empty, pGuard = Map.empty, pGuardFunc = Map.empty, pTile = Map.empty, pTileOrder = 0, pColorNum = 0 }

-- icc = "g++"
icc = "icpc"
-- icc = "g++-cilk"

-- compilation flags for icpc
iccFlags = ["-O3", "-DNDEBUG", "-std=c++0x", "-Wall", "-Werror", "-ipo"]
-- compilation flags for g++
-- iccFlags = ["-O3", "-DNDEBUG", "-std=c++0x", "-Wall", "-Werror"]

-- preprocessing flags for g++
-- iccPPFlags = ["-E", "-DNCHECK_SHAPE", "-DNDEBUG", "-std=c++0x", "-Wall", "-Werror"]
-- preprocessing flags for icpc
iccPPFlags = ["-P", "-C", "-DNCHECK_SHAPE", "-DNDEBUG", "-std=c++0x", "-Wall", "-Werror", "-ipo"]

-- iccDebugPPFlags = ["-P", "-C", "-DCHECK_SHAPE", "-DDEBUG", "-g3", "-std=c++0x", "-include", "cilk_stub.h"]
-- iccDebugPPFlags = ["-P", "-C", "-DCHECK_SHAPE", "-DDEBUG", "-g3", "-std=c++0x"]
-- debug flags for icpc
iccDebugPPFlags = ["-P", "-C", "-DCHECK_SHAPE", "-DNDEBUG", "-g3", "-std=c++0x"]
-- debug flags for g++
-- iccDebugPPFlags = ["-E", "-DCHECK_SHAPE", "-DNDEBUG", "-g3", "-std=c++0x"]

parseArgs :: ([String], [String], PMode, Bool, Bool, [String]) -> [String] -> ([String], [String], PMode, Bool, Bool, [String])
parseArgs (inFiles, inDirs, mode, debug, showFile, userArgs) aL 
    | elem "--help" aL =
        let l_mode = PHelp
            aL' = delete "--help" aL
        in  (inFiles, inDirs, l_mode, debug, showFile, aL')
    | elem "-h" aL =
        let l_mode = PHelp
            aL' = delete "-h" aL
        in  (inFiles, inDirs, l_mode, debug, showFile, aL')
    | elem "-auto-optimize" aL =
        let l_mode = PDefault
            aL' = delete "-auto-optimize" aL
        in  (inFiles, inDirs, l_mode, debug, showFile, aL')
------------------------------------------------------------------------------
-- so far, the split-caching mode doesn't work for multiple-kernel case!!!! --
------------------------------------------------------------------------------
    | elem "-split-caching" aL =
        let l_mode = PCaching
            aL' = delete "-split-caching" aL
        in  parseArgs (inFiles, inDirs, l_mode, debug, showFile, aL') aL'
    | elem "-all-cond-tile-macro-overlap" aL =
        let l_mode = PAllCondTileMacroOverlap
            aL' = delete "-all-cond-tile-macro-overlap" aL
        in  parseArgs (inFiles, inDirs, l_mode, debug, showFile, aL') aL'
    | elem "-all-cond-tile-c-pointer-overlap" aL =
        let l_mode = PAllCondTileCPointerOverlap
            aL' = delete "-all-cond-tile-c-pointer-overlap" aL
        in  parseArgs (inFiles, inDirs, l_mode, debug, showFile, aL') aL'
    | elem "-all-cond-tile-pointer-overlap" aL =
        let l_mode = PAllCondTilePointerOverlap
            aL' = delete "-all-cond-tile-pointer-overlap" aL
        in  parseArgs (inFiles, inDirs, l_mode, debug, showFile, aL') aL'
    | elem "-all-cond-tile-opt-pointer-overlap" aL =
        let l_mode = PAllCondTileOptPointerOverlap
            aL' = delete "-all-cond-tile-opt-pointer-overlap" aL
        in  parseArgs (inFiles, inDirs, l_mode, debug, showFile, aL') aL'
    | elem "-unroll-t-tile-macro-overlap" aL =
        let l_mode = PUnrollTimeTileMacroOverlap
            aL' = delete "-unroll-t-tile-macro-overlap" aL
        in  parseArgs (inFiles, inDirs, l_mode, debug, showFile, aL') aL'
    | elem "-unroll-t-tile-c-pointer-overlap" aL =
        let l_mode = PUnrollTimeTileCPointerOverlap
            aL' = delete "-unroll-t-tile-c-pointer-overlap" aL
        in  parseArgs (inFiles, inDirs, l_mode, debug, showFile, aL') aL'
    | elem "-unroll-t-tile-pointer-overlap" aL =
        let l_mode = PUnrollTimeTilePointerOverlap
            aL' = delete "-unroll-t-tile-pointer-overlap" aL
        in  parseArgs (inFiles, inDirs, l_mode, debug, showFile, aL') aL'
    | elem "-unroll-t-tile-opt-pointer-overlap" aL =
        let l_mode = PUnrollTimeTileOptPointerOverlap
            aL' = delete "-unroll-t-tile-opt-pointer-overlap" aL
        in  parseArgs (inFiles, inDirs, l_mode, debug, showFile, aL') aL'
    | elem "-all-cond-tile-macro" aL =
        let l_mode = PAllCondTileMacro
            aL' = delete "-all-cond-tile-macro" aL
        in  parseArgs (inFiles, inDirs, l_mode, debug, showFile, aL') aL'
    | elem "-all-cond-tile-c-pointer" aL =
        let l_mode = PAllCondTileCPointer
            aL' = delete "-all-cond-tile-c-pointer" aL
        in  parseArgs (inFiles, inDirs, l_mode, debug, showFile, aL') aL'
    | elem "-all-cond-tile-pointer" aL =
        let l_mode = PAllCondTilePointer
            aL' = delete "-all-cond-tile-pointer" aL
        in  parseArgs (inFiles, inDirs, l_mode, debug, showFile, aL') aL'
    | elem "-all-cond-tile-opt-pointer" aL =
        let l_mode = PAllCondTileOptPointer
            aL' = delete "-all-cond-tile-opt-pointer" aL
        in  parseArgs (inFiles, inDirs, l_mode, debug, showFile, aL') aL'
    | elem "-unroll-t-tile-macro" aL =
        let l_mode = PUnrollTimeTileMacro
            aL' = delete "-unroll-t-tile-macro" aL
        in  parseArgs (inFiles, inDirs, l_mode, debug, showFile, aL') aL'
    | elem "-unroll-t-tile-c-pointer" aL =
        let l_mode = PUnrollTimeTileCPointer
            aL' = delete "-unroll-t-tile-c-pointer" aL
        in  parseArgs (inFiles, inDirs, l_mode, debug, showFile, aL') aL'
    | elem "-unroll-t-tile-pointer" aL =
        let l_mode = PUnrollTimeTilePointer
            aL' = delete "-unroll-t-tile-pointer" aL
        in  parseArgs (inFiles, inDirs, l_mode, debug, showFile, aL') aL'
    | elem "-unroll-t-tile-opt-pointer" aL =
        let l_mode = PUnrollTimeTileOptPointer
            aL' = delete "-unroll-t-tile-opt-pointer" aL
        in  parseArgs (inFiles, inDirs, l_mode, debug, showFile, aL') aL'
    | elem "-unroll-multi-kernel" aL =
        let l_mode = PMUnroll
            aL' = delete "-unroll-multi-kernel" aL
        in  parseArgs (inFiles, inDirs, l_mode, debug, showFile, aL') aL'
    | elem "-split-c-pointer" aL =
        let l_mode = PCPointer
            aL' = delete "-split-c-pointer" aL
        in  parseArgs (inFiles, inDirs, l_mode, debug, showFile, aL') aL'
    | elem "-split-opt-pointer" aL =
        let l_mode = POptPointer
            aL' = delete "-split-opt-pointer" aL
        in  parseArgs (inFiles, inDirs, l_mode, debug, showFile, aL') aL'
    | elem "-split-pointer" aL =
        let l_mode = PPointer
            aL' = delete "-split-pointer" aL
        in  parseArgs (inFiles, inDirs, l_mode, debug, showFile, aL') aL'
    | elem "-split-macro-shadow" aL =
        let l_mode = PMacroShadow
            aL' = delete "-split-macro-shadow" aL
        in  parseArgs (inFiles, inDirs, l_mode, debug, showFile, aL') aL'
    | elem "-showFile" aL =
        let l_showFile = True
            aL' = delete "-showFile" aL
        in  parseArgs (inFiles, inDirs, mode, debug, l_showFile, aL') aL'
    | elem "-debug" aL =
        let l_debug = True 
            l_mode = PDebug
            aL' = delete "-debug" aL
        in  parseArgs (inFiles, inDirs, l_mode, l_debug, showFile, aL') aL'
    | null aL == False =
        let (l_files, l_dirs, l_mode, aL') = findCPP aL ([], [], mode, aL)
        in  (l_files, l_dirs, l_mode, debug, showFile, aL')
    | otherwise = 
        let l_mode = PNoPP
        in  (inFiles, inDirs, l_mode, debug, showFile, aL)

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


