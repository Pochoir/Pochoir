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

-- We don't compact pRegLambdaTerm because of the l_run_kernel
pRegLambdaTerm :: PMode -> Int -> PStencil -> (PGuard, [PTile]) -> String
pRegLambdaTerm l_mode l_rank l_stencil (g, t) =
    let l_tag = getTagFromMode l_mode
        l_regBound = sRegBound l_stencil
        l_order = gOrder g
        l_bdry_name = l_tag ++ "_boundary_kernel_" ++ show l_order
        l_cond_bdry_name = l_tag ++ "_cond_boundary_kernel_" ++ show l_order
        l_obase_name = l_tag ++ "_interior_kernel_" ++ show l_order
        l_cond_obase_name = l_tag ++ "_cond_interior_kernel_" ++ show l_order
        l_bdry_pointer = (mkParen . derefPointer) l_bdry_name
        l_cond_bdry_pointer = (mkParen . derefPointer) l_cond_bdry_name
        l_obase_pointer = (mkParen . derefPointer) l_obase_name
        l_cond_obase_pointer = (mkParen . derefPointer) l_cond_obase_name
        l_run_kernel = 
            if l_mode == PAllCondTileMacroOverlap ||
               l_mode == PAllCondTileCPointerOverlap ||
               l_mode == PAllCondTilePointerOverlap ||
               l_mode == PAllCondTileOptPointerOverlap 
               then if l_regBound 
                       then l_obase_pointer ++ ", " ++ l_bdry_pointer
                       else l_obase_pointer
               else if l_regBound
                       then l_obase_pointer ++ ", " ++ l_cond_obase_pointer ++
                            ", " ++ l_bdry_pointer ++ ", " ++ l_cond_bdry_pointer
                       else l_obase_pointer ++ ", " ++ l_cond_obase_pointer
        l_pochoir_id = sName l_stencil
        l_guard_name = pGetOverlapGuardName g
        l_unroll = foldr max 0 $ map pTileLength t
    in  pTab ++ l_pochoir_id ++ ".Register_Tile_Obase_Kernels(" ++ l_guard_name ++
        ", " ++ show l_unroll ++ ", " ++ l_run_kernel ++ ");" ++ breakline

pDestroyLambdaTerm :: PMode -> Int -> String -> PStencil -> (PGuard, [PTile]) -> String
pDestroyLambdaTerm l_mode l_rank l_str l_stencil (g, t) =
    let l_tag = getTagFromMode l_mode
        l_regBound = sRegBound l_stencil
        l_order = gOrder g
        l_name = l_tag ++ "_" ++ l_str ++ "_kernel_" ++ show l_order
        l_pointer = mkInput l_name
        l_del_lambdaPointer =
            if (l_regBound && l_str == "boundary") || 
               (l_regBound && l_str == "cond_boundary") ||
               (l_str == "interior") ||
               (l_str == "cond_interior")
               then pDelPointer  l_pointer 
               else ""
        l_del_kernelPointer =
            if (l_regBound && l_str == "boundary") || 
               (l_regBound && l_str == "cond_boundary") ||
               (l_str == "interior") ||
               (l_str == "cond_interior")
               then pDelPointer  l_name
               else ""
    in  breakline ++ l_del_lambdaPointer ++
        breakline ++ l_del_kernelPointer 

pCreateLambdaTerm :: PMode -> Int -> String -> PStencil -> [String] -> (PGuard, [PTile]) -> String
pCreateLambdaTerm l_mode l_rank l_str l_stencil l_inputParams (g, t) =
    let l_tag = getTagFromMode l_mode
        l_regBound = sRegBound l_stencil
        l_order = gOrder g
        l_shape_name = "__POCHOIR_Shape__" ++ show l_order ++ "__"
        l_name = l_tag ++ "_" ++ l_str ++ "_kernel_" ++ show l_order
        l_class = pSys l_name
        l_pointer = mkInput l_name
        l_kernel_pointer = mkLocal l_name
        l_kernel_class = "Pochoir_Obase_Kernel < " ++ 
                         "decltype" ++ (mkParen $ derefPointer l_pointer) ++
                         ", " ++ show l_rank ++ " > "
        l_kernel_inputs = [(mkParen . derefPointer) l_pointer] ++ [l_shape_name] 
        l_new_lambdaPointer =
            if (l_regBound && l_str == "boundary") || 
               (l_regBound && l_str == "cond_boundary") ||
               (l_str == "interior") ||
               (l_str == "cond_interior")
               then pNewPointer  l_pointer l_class l_inputParams
               else "" 
        l_new_kernelPointer =
            if (l_regBound && l_str == "boundary") || 
               (l_regBound && l_str == "cond_boundary") ||
               (l_str == "interior") ||
               (l_str == "cond_interior")
               then pNewPointer  l_kernel_pointer l_kernel_class l_kernel_inputs
               else "" 
        l_decl_local_obase = pDeclVar (mkPointer l_kernel_class) l_kernel_pointer "NULL"
        l_assign = pAssign l_name l_kernel_pointer
    in  breakline ++ pTab ++ l_new_lambdaPointer ++
        breakline ++ pTab ++ l_decl_local_obase ++
        breakline ++ pTab ++ l_new_kernelPointer ++
        breakline ++ pTab ++ l_assign ++ breakline

{-
 - a version using C-style calloc()/free()
pNewPointer :: String -> String -> [String] -> String
pNewPointer l_pointer l_class l_inputParams =
    let l_alloc = l_pointer ++ " = " ++ mkPointer l_class ++ "calloc(1, " ++ 
                  "sizeof" ++ mkParen l_class ++ ");"
        l_init_content = 
                breakline ++ "if ( " ++ l_pointer ++ " == NULL ) {" ++
                breakline ++ pTab ++ "printf(\" Failure in create_lambda allocation!\\n\");" ++
                breakline ++ pTab ++ "exit(EXIT_FAILURE);" ++ breakline ++ "}" ++ 
                l_pointer ++ "->Init" ++ (mkParen $ intercalate ", " l_inputParams) ++ ";"
    in  l_alloc ++ breakline ++ l_init_content
 -}
 -- a version using C++ style new()/delete()
pNewPointer :: String -> String -> [String] -> String
pNewPointer l_pointer l_class l_inputParams =
    let l_alloc = l_pointer ++ " = new " ++ l_class ++ 
                  (mkParen $ intercalate ", " l_inputParams) ++ ";" 
    in  l_alloc 

{-
 - a version using C-style calloc()/free()
pDelPointer :: String -> String
pDelPointer l_pointer =
    let l_dealloc = 
            breakline ++ "if " ++ mkParen (l_pointer ++ " != NULL") ++ " {" ++
            breakline ++ pTab ++ "free " ++ mkParen l_pointer ++ ";" ++ breakline ++ "}"
    in  l_dealloc 
 -}
 -- a version using C++ style new()/delete()
pDelPointer :: String -> String
pDelPointer l_pointer =
    let l_dealloc = pTab ++ "delete " ++ l_pointer ++ ";"
    in  l_dealloc 
