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
import Data.Maybe
import System.Directory 
import System.Cmd (rawSystem)
import Data.Char (isSpace)
import qualified Data.Map as Map
import Text.ParserCombinators.Parsec 

import PGenData
import PGenMainParser
import PGenBasicParser
import PGenUtils

main :: IO ()
main = do args <- getArgs
          whilst (null args || length args < 2) $ do
             printUsage
             exitFailure
          let (stencil_name, color_file, kernel_file, mode, showFile, userArgs) 
                = parseArgs ("", "", "", PDefault, True, []) args
          whilst (mode == PHelp) $ do
             printOptions
             exitFailure
          whilst (mode /= PNoPP) $ do
             putStrLn ("stencil_name = " ++ stencil_name ++
                       ", color_file = " ++ color_file ++
                       ", kernel_file = " ++ kernel_file ++ ", mode = " ++ show mode)
             colorh <- openFile color_file ReadMode
             colorVectors <- pParseColorFile colorh
             hClose colorh
             kernelh <- openFile kernel_file ReadMode
             -- let genKernel_file = pSubstitute "kernel_info" "gen_kernel" kernel_file
             let genKernel_file = pSubstitute "color.dat" "gen_kernel.cpp" color_file
             genKernelh <- openFile genKernel_file WriteMode
             pGenKernel mode stencil_name colorVectors kernelh genKernelh
             hClose genKernelh
          whilst (showFile == False) $ do
             removeFile color_file 
             removeFile kernel_file

whilst :: Bool -> IO () -> IO ()
whilst True action = action
whilst False action = return () 

pInitState = ParserState { pMode = PUnrollPointer, pColorVectors = [], pArray = Map.empty, pStencil = Map.empty, pShape = Map.empty, pRange = Map.empty, pKernel = Map.empty, pKernelFunc = Map.empty, pGuard = Map.empty, pGuardFunc = Map.empty, pTile = Map.empty, pTileOrder = 0, pGuardOrder = 0, pGuardFuncOrder = 0, pStencilName = ""}

-- icc = "icpc"
-- icc = "g++_4_7"

-- iccFlags = ["-O3", "-DNDEBUG", "-std=c++11", "-Wall", "-Werror", "-ipo", "-fcilkplus", "-lcilkrts"]

-- iccPPFlags = ["-E", "-C", "-DNCHECK_SHAPE", "-DNDEBUG", "-std=c++11", "-Wall", "-Werror", "-ipo", "-fcilkplus"]

-- iccDebugFlags = ["-DDEBUG", "-O0", "-g3", "-std=c++0x", "-include", "cilk_stub.h"]
-- iccDebugFlags = ["-DNDEBUG", "-O0", "-g3", "-std=c++11", "-include", "cilk_stub.h"]

-- iccDebugPPFlags = ["-P", "-C", "-DCHECK_SHAPE", "-DDEBUG", "-g3", "-std=c++0x", "-include", "cilk_stub.h"]
-- iccDebugPPFlags = ["-P", "-C", "-DCHECK_SHAPE", "-DDEBUG", "-g3", "-std=c++0x"]
-- iccDebugPPFlags = ["-E", "-C", "-DNDEBUG", "-g3", "-std=c++11", "-include", "cilk_stub.h"]

parseArgs :: (String, String, String, PMode, Bool, [String]) -> [String] -> (String, String, String, PMode, Bool, [String])
parseArgs (stencil_name, color_file, kernel_file, mode, showFile, userArgs) aL 
    | elem "--help" aL =
        let l_mode = PHelp
            aL' = delete "--help" aL
        in  (stencil_name, color_file, kernel_file, l_mode, showFile, aL')
    | elem "-h" aL =
        let l_mode = PHelp
            aL' = delete "-h" aL
        in  (stencil_name, color_file, kernel_file, l_mode, showFile, aL')
    | elem "-showFile" aL =
        let l_showFile = True
            aL' = delete "-showFile" aL
        in  parseArgs (stencil_name, color_file, kernel_file, mode, l_showFile, aL') aL'
    | elem "-auto-optimize" aL =
        let l_mode = PDefault
            aL' = delete "-auto-optimize" aL
        in  parseArgs (stencil_name, color_file, kernel_file, l_mode, showFile, aL') aL'
------------------------------------------------------------------------------
-- so far, the split-caching mode doesn't work for multiple-kernel case!!!! --
------------------------------------------------------------------------------
    | elem "-split-caching" aL =
        let l_mode = PCaching
            aL' = delete "-split-caching" aL
        in  parseArgs (stencil_name, color_file, kernel_file, l_mode, showFile, aL') aL'
    | elem "-cond-macro" aL =
        let l_mode = PCondMacro
            aL' = delete "-cond-macro" aL
        in  parseArgs (stencil_name, color_file, kernel_file, l_mode, showFile, aL') aL'
    | elem "-cond-c-pointer" aL =
        let l_mode = PCondCPointer
            aL' = delete "-cond-c-pointer" aL
        in  parseArgs (stencil_name, color_file, kernel_file, l_mode, showFile, aL') aL'
    | elem "-cond-pointer" aL =
        let l_mode = PCondPointer
            aL' = delete "-cond-pointer" aL
        in  parseArgs (stencil_name, color_file, kernel_file, l_mode, showFile, aL') aL'
    | elem "-cond-opt-pointer" aL =
        let l_mode = PCondOptPointer
            aL' = delete "-cond-opt-pointer" aL
        in  parseArgs (stencil_name, color_file, kernel_file, l_mode, showFile, aL') aL'
    | elem "-unroll-macro" aL =
        let l_mode = PUnrollMacro
            aL' = delete "-unroll-macro" aL
        in  parseArgs (stencil_name, color_file, kernel_file, l_mode, showFile, aL') aL'
    | elem "-unroll-c-pointer" aL =
        let l_mode = PUnrollCPointer
            aL' = delete "-unroll-c-pointer" aL
        in  parseArgs (stencil_name, color_file, kernel_file, l_mode, showFile, aL') aL'
    | elem "-unroll-pointer" aL =
        let l_mode = PUnrollPointer
            aL' = delete "-unroll-pointer" aL
        in  parseArgs (stencil_name, color_file, kernel_file, l_mode, showFile, aL') aL'
    | elem "-unroll-opt-pointer" aL =
        let l_mode = PUnrollOptPointer
            aL' = delete "-unroll-opt-pointer" aL
        in  parseArgs (stencil_name, color_file, kernel_file, l_mode, showFile, aL') aL'
    | elem "-unroll-multi-kernel" aL =
        let l_mode = PMUnroll
            aL' = delete "-unroll-multi-kernel" aL
        in  parseArgs (stencil_name, color_file, kernel_file, l_mode, showFile, aL') aL'
    | elem "-stencil" aL = 
        let idx = fromMaybe 0 $ elemIndex "-stencil" aL
            l_stencil_name = aL !! (idx + 1)
            aL' = delete "-stencil" aL
            aL'' = delete l_stencil_name aL'
        in  parseArgs (l_stencil_name, color_file, kernel_file, mode, showFile, aL'') aL''
    | null aL == False =
        let (l_color_file, l_color_dir) = findFileBySuffix aL "color.dat" 
            (l_kernel_file, l_kernel_dir) = findFileBySuffix aL "kernel_info.cpp" 
            aL' = delete (l_color_dir ++ l_color_file) $ delete (l_kernel_dir ++ l_kernel_file) aL
        in  (stencil_name, l_color_dir++l_color_file, l_kernel_dir++l_kernel_file, mode, showFile, aL')
    | otherwise = 
        let l_mode = PNoPP
        in  (stencil_name, color_file, kernel_file, l_mode, showFile, aL)

findFileBySuffix :: [String] -> String -> (String, String) 
findFileBySuffix [] _ = ("", "")
findFileBySuffix (a:as) suffix
    | isSuffixOf suffix a = 
        let l_file = drop (1 + (pLast $ findIndices (== '/') a)) a
            l_dir  = take (1 + (pLast $ findIndices (== '/') a)) a
            pLast [] = -1
            pLast aL@(a:as) = last aL
            -- l_pochoir_file = rename "_pochoir" l_file 
        in  (l_file, l_dir)
    | otherwise = findFileBySuffix as suffix

printUsage :: IO ()
printUsage =
    do putStrLn ("Usage: genkernels [OPTION] [KERNEL_INFO_FILE] [COLOR_INFO_FILE]")
       putStrLn ("Try `genkernels --help' for more options.")

printOptions :: IO ()
printOptions = 
    do putStrLn ("Usage: genkernels [OPTION] [KERNEL_INFO_FILE] [COLOR_INFO_FILE]")
       putStrLn ("Run the genkernels with kernel_info and color_info.")

pGenKernel :: PMode -> String -> [Homogeneity] -> Handle -> Handle -> IO ()           
pGenKernel mode stencil_name colorVectors kernelh genKernelh = 
    do ls <- hGetContents kernelh
       let pRevInitState = pInitState { pMode = mode, pColorVectors = colorVectors, pStencilName = stencil_name}
       case runParser pParser pRevInitState "" $ stripWhite ls of
           Left err -> print err
           Right str -> hPutStrLn genKernelh str

pParseColorFile :: Handle -> IO [Homogeneity]
pParseColorFile colorh =
    do ls <- hGetContents colorh
       case parse pParseColor "" ls of
           Left err -> 
                do print err
                   return []
           Right colorVectors -> 
                do putStrLn $ show colorVectors
                   return colorVectors
           
pParseColor :: Parser [Homogeneity]
pParseColor = many $ parens pColor

pColor :: Parser Homogeneity
pColor = 
    do orStr <- many digit 
       comma
       andStr <- many digit 
       let l_homo = Homogeneity { size = length orStr, o = pBinToDec orStr, a = pBinToDec andStr }
       return l_homo

