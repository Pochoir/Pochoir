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

-- The Basic Parser for a second pass --
module PGenBasicParser2 where

import Text.ParserCombinators.Parsec
import qualified Text.ParserCombinators.Parsec.Token as Token
import Text.ParserCombinators.Parsec.Expr
import Text.ParserCombinators.Parsec.Language
import Text.Read (read)

import Data.Char
import Data.List
import qualified Data.Map as Map

import PGenBasicParser
import PGenShow
import PGenUtils
import PGenData

pRegLambdaTerm :: PMode -> Int -> PStencil -> (PGuard, [PTile]) -> String
pRegLambdaTerm l_mode l_rank l_stencil (g, t) =
    let l_tag = getTagFromMode l_mode
        l_regBound = sRegBound l_stencil
        l_order = gOrder g
        l_bdry_name = l_tag ++ "_boundary_kernel_" ++ show l_order
        l_obase_name = l_tag ++ "_interior_kernel_" ++ show l_order
        l_bdry_pointer = (mkParen . derefPointer) l_bdry_name
        l_obase_pointer = (mkParen . derefPointer) l_obase_name
        l_run_kernel = if l_regBound 
                          then l_obase_pointer ++ ", " ++ l_bdry_pointer
                          else l_obase_pointer
        l_pochoir_id = sName l_stencil
        l_guard_name = pGetOverlapGuardName g
        l_unroll = foldr max 0 $ map pTileLength t
    in  l_pochoir_id ++ ".Register_Tile_Obase_Kernels(" ++ l_guard_name ++
        ", " ++ show l_unroll ++ ", " ++ l_run_kernel ++ ");" ++ breakline

pDestroyLambdaTerm :: PMode -> Int -> PStencil -> (PGuard, [PTile]) -> String
pDestroyLambdaTerm l_mode l_rank l_stencil (g, t) =
    let l_tag = getTagFromMode l_mode
        l_regBound = sRegBound l_stencil
        l_order = gOrder g
        l_bdry_name = l_tag ++ "_boundary_kernel_" ++ show l_order
        l_obase_name = l_tag ++ "_interior_kernel_" ++ show l_order
        l_bdry_pointer = mkInput l_bdry_name
        l_obase_pointer = mkInput l_obase_name
        l_del_bdry_lambdaPointer =
            if l_regBound
               then pDelPointer l_bdry_pointer 
               else ""
        l_del_obase_lambdaPointer = pDelPointer l_obase_pointer 
        l_del_bdry_kernelPointer =
            if l_regBound
               then pDelPointer l_bdry_name
               else ""
        l_del_obase_kernelPointer = pDelPointer l_obase_name
    in  breakline ++ l_del_bdry_lambdaPointer ++
        breakline ++ l_del_bdry_kernelPointer ++
        breakline ++ l_del_obase_lambdaPointer ++
        breakline ++ l_del_obase_kernelPointer

pCreateLambdaTerm :: PMode -> Int -> PStencil -> [String] -> (PGuard, [PTile]) -> String
pCreateLambdaTerm l_mode l_rank l_stencil l_inputParams (g, t) =
    let l_tag = getTagFromMode l_mode
        l_regBound = sRegBound l_stencil
        l_order = gOrder g
        l_shape_name = "__POCHOIR_Shape__" ++ show l_order ++ "__"
        l_bdry_name = l_tag ++ "_boundary_kernel_" ++ show l_order
        l_obase_name = l_tag ++ "_interior_kernel_" ++ show l_order
        l_bdry_class = pSys l_bdry_name
        l_obase_class = pSys l_obase_name
        l_bdry_pointer = mkInput l_bdry_name
        l_obase_pointer = mkInput l_obase_name
        l_bdry_kernel_pointer = l_bdry_name
        l_obase_kernel_pointer = l_obase_name
        l_kernel_class = "Pochoir_Obase_Kernel < " ++ show l_rank ++ " > "
        l_bdry_kernel_inputs = [(mkParen . derefPointer) l_bdry_pointer] ++ [l_shape_name] 
        l_obase_kernel_inputs = [(mkParen . derefPointer) l_obase_pointer] ++ [l_shape_name] 
        l_new_bdry_lambdaPointer =
            if l_regBound
               then pNewLambdaClosure l_bdry_pointer l_bdry_class l_inputParams
               else ""
        l_new_obase_lambdaPointer = pNewLambdaClosure l_obase_pointer l_obase_class l_inputParams
        l_new_bdry_kernelPointer =
            if l_regBound
               then pNewLambdaClosure l_bdry_kernel_pointer l_kernel_class l_bdry_kernel_inputs
               else ""
        l_new_obase_kernelPointer = pNewLambdaClosure l_obase_kernel_pointer l_kernel_class l_obase_kernel_inputs
    in  breakline ++ l_new_bdry_lambdaPointer ++
        breakline ++ l_new_bdry_kernelPointer ++
        breakline ++ l_new_obase_lambdaPointer ++
        breakline ++ l_new_obase_kernelPointer

pNewLambdaClosure :: String -> String -> [String] -> String
pNewLambdaClosure l_pointer l_class l_inputParams =
    l_pointer ++ " = new (std::nothrow)" ++ l_class ++
    -- l_pointer ++ " = " ++ mkPointer l_class ++ "calloc(1, " ++ l_class ++ ")" ++
    (mkParen $ intercalate ", " l_inputParams) ++ ";" ++
    breakline ++ "if ( " ++ l_pointer ++
    " == NULL ) {" ++
    breakline ++ pTab ++ "printf(\" Failure in create_lambda allocation!\\n\");" ++
    breakline ++ pTab ++ "exit(EXIT_FAILURE);" ++ breakline ++ "}"

pDelPointer :: String -> String
pDelPointer l_pointer =
    breakline ++ "if " ++ mkParen (l_pointer ++ " != NULL") ++ " {" ++ 
    breakline ++ pTab ++ "delete " ++ l_pointer ++ ";" ++ 
    -- breakline ++ pTab ++ "free " ++ mkParen l_pointer ++ ";" ++ 
    breakline ++ pTab ++ l_pointer ++ " = NULL;" ++
    breakline ++ "}"
