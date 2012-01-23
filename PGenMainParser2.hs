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

-- The Main Parser for a second pass --
module PGenMainParser2 where

import Text.ParserCombinators.Parsec

import Control.Monad

import PGenBasicParser
import PGenBasicParser2
import PGenUtils
import PGenData
import PGenShow
-- import Text.Show
import Data.List
import Data.Bits
import qualified Data.Map as Map

pCodeGen :: PMode -> [Homogeneity] -> [PGuardFunc] -> PStencil -> (String, [(PGuard, [PTile])])
pCodeGen l_mode l_color_vectors l_guardFuncs l_stencil = 
    let -- we are assuming that all guards and kernels are of the same rank
        l_reg_GTs = sRegTileKernel l_stencil
        l_rank = if length l_reg_GTs > 0
                    then gRank $ fst $ head l_reg_GTs
                    else 0
        l_guardTiles = pGetGuardTiles l_mode 0 l_rank l_color_vectors l_reg_GTs
        l_guards = map (pShowGlobalGuardString " && ") l_guardTiles
        l_tiles = map (pShowAutoTileString l_mode l_stencil) l_guardTiles
        l_create_lambdas = pCreateLambdas l_mode l_rank l_stencil l_guardTiles
        l_destroy_lambdas = pDestroyLambdas l_mode l_rank l_stencil l_guardTiles
        l_reg_lambdas = pRegLambdas l_mode l_rank l_stencil l_guardTiles
        l_str_GTs = zipWith (++) l_guards l_tiles
    in  (breakline ++ pShowHeader ++ 
            breakline ++ "/* Original Codes */" ++
            breakline ++ pShowGuardFuncList l_guardFuncs ++
            breakline ++ "/* Generated Codes */" ++
            breakline ++ pShowColorVectors l_color_vectors ++ 
            breakline ++ concat l_str_GTs ++
            breakline ++ l_create_lambdas ++
            breakline ++ l_destroy_lambdas ++
            breakline ++ l_reg_lambdas
            , l_guardTiles)

pRegLambdas :: PMode -> Int -> PStencil -> [(PGuard, [PTile])] -> String
pRegLambdas _ _ _ [] = ""
pRegLambdas l_mode l_rank l_stencil cL =
    let l_pochoir_id = sName l_stencil
        l_pochoir_ref = "Pochoir < " ++ show l_rank ++ " > & "
        l_header = externC ++ "int Register_Lambdas " ++
                   mkParen (l_pochoir_ref ++ l_pochoir_id) ++ " { " ++ breakline
        l_content = concatMap (pRegLambdaTerm l_mode l_rank l_stencil) cL
        l_tail = breakline ++ "return 0;" ++ breakline ++ "}"
    in  l_header ++ l_content ++ l_tail

pDestroyLambdas :: PMode -> Int -> PStencil -> [(PGuard, [PTile])] -> String
pDestroyLambdas _ _ _ [] = ""
pDestroyLambdas l_mode l_rank l_stencil cL =
    let l_header = externC ++ "int Destroy_Lambdas " ++
                   mkParen "void" ++ " { " ++ breakline
        l_content = concatMap (pDestroyLambdaTerm l_mode l_rank l_stencil) cL
        l_tail = breakline ++ "return 0;" ++ breakline ++ "}"
    in  l_header ++ l_content ++ l_tail

pCreateLambdas :: PMode -> Int -> PStencil -> [(PGuard, [PTile])] -> String
pCreateLambdas _ _ _ [] = ""
pCreateLambdas l_mode l_rank l_stencil cL =
    let l_arrayInUse = sArrayInUse l_stencil
        l_arrayList = map aName l_arrayInUse
        l_arrayInputList = map (mkInput . aName) l_arrayInUse
        l_arrayInputRefList = map pShowPochoirArrayRef $ zip3 (map aRank l_arrayInUse) (map aType l_arrayInUse) (map aName l_arrayInUse)
        l_header = externC ++ "int Create_Lambdas " ++ 
                    mkParen (intercalate ", " l_arrayInputRefList) ++ "{" ++ breakline
        l_content = concatMap (pCreateLambdaTerm l_mode l_rank l_stencil l_arrayInputList) cL
        l_tail = "}\n"
    in  l_header ++ l_content ++ l_tail

pShowGuardFuncList :: [PGuardFunc] -> String
pShowGuardFuncList [] = ""
pShowGuardFuncList gL@(g:gs) = pShowGlobalGuardFunc (gfName g) g ++ pShowGuardFuncList gs

pShowColorVectors :: [Homogeneity] -> String 
pShowColorVectors l_color_vectors = 
    breakline ++ "/* " ++ 
    breakline ++ show l_color_vectors ++ 
    breakline ++ " */" ++ breakline

pGetGuardTiles :: PMode -> Int -> Int -> [Homogeneity] -> [(PGuard, PTile)] -> [(PGuard, [PTile])]
pGetGuardTiles _ _ _ [] _ = []
pGetGuardTiles l_mode l_order l_rank cL@(c:cs) l_reg_GTs =
    let l_gts = pGetGuardTilesTerm l_mode l_order l_rank c l_reg_GTs  
    in  [l_gts] ++ pGetGuardTiles l_mode (l_order+1) l_rank cs l_reg_GTs

pGetGuardTilesTerm :: PMode -> Int -> Int -> Homogeneity -> [(PGuard, PTile)] -> (PGuard, [PTile])
pGetGuardTilesTerm l_mode l_order l_rank l_color l_reg_GTs =
    let l_guards = map fst l_reg_GTs
        l_guards' = zipWith pFillGuardOrder [0..] l_guards
        l_out_guard = foldr (pColorGuard l_rank l_color) emptyGuard l_guards'
        l_out_guard' = pFillGuardOrder l_order l_out_guard
        l_tiles = map snd l_reg_GTs
        l_tiles' = zipWith pFillTileOrder [0..] l_tiles
        l_out_tiles = foldr (pColorTile l_rank l_color) [] l_tiles'
        l_out_tiles' = map (pFillTileOrder l_order) l_out_tiles
    in  (l_out_guard', l_out_tiles')

pColorGuard :: Int -> Homogeneity -> PGuard -> PGuard -> PGuard
pColorGuard l_rank l_color l_old_guard_l l_old_guard_r =
    let l_idx = size l_color - 1 - gOrder l_old_guard_l
        l_condStr = 
            if testBit (o l_color) l_idx && testBit (a l_color) l_idx
               then [gName l_old_guard_l] ++ gComment l_old_guard_r 
               else if (not $ testBit (o l_color) l_idx) && 
                        (not $ testBit (a l_color) l_idx)
                       then ["!" ++ gName l_old_guard_l] ++ gComment l_old_guard_r
                       else gComment l_old_guard_r
    in  PGuard { gName = pSys $ intercalate "_" l_condStr, gRank = l_rank, gFunc = emptyGuardFunc, gComment = l_condStr, gOrder = 0, gColor = l_color }
                            
pColorTile :: Int -> Homogeneity -> PTile -> [PTile] -> [PTile]
pColorTile l_rank l_color l_old_tile_l l_old_tiles_r =
    let l_idx = size l_color - 1 - tOrder l_old_tile_l
        l_tiles = 
            if testBit (o l_color) l_idx && testBit (a l_color) l_idx
               then let l_new_tile_l = pSetTileOp PSERIAL l_old_tile_l
                        l_new_tile_l' = l_new_tile_l { tColor = l_color }
                    in  [l_new_tile_l'] ++ l_old_tiles_r
               else if (testBit (o l_color) l_idx) &&
                        (not $ testBit (a l_color) l_idx)
                       then let l_new_tile_l = pSetTileOp PINCLUSIVE l_old_tile_l
                                l_new_tile_l' = l_new_tile_l { tColor = l_color }
                            in  [l_new_tile_l'] ++ l_old_tiles_r
                       else l_old_tiles_r
    in  l_tiles

pShowHeader :: String
pShowHeader = 
    breakline ++ "#include <cstdio>" ++
    breakline ++ "#include <cstdlib>" ++
    breakline ++ "#include <cassert>" ++
    breakline ++ "#include <functional>" ++
    breakline ++ "#include <dlfcn.h>" ++ 
    breakline 
