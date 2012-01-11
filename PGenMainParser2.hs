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

pCodeGen :: PMode -> [Homogeneity] -> PStencil -> (String, [(PGuard, [PTile])])
pCodeGen l_mode l_color_vectors l_stencil = 
    let -- we are assuming that all guards and kernels are of the same rank
        l_reg_GTs = sRegTileKernel l_stencil
        l_rank = if length l_reg_GTs > 0
                    then gRank $ fst $ head l_reg_GTs
                    else 0
        l_guardTiles = pGetGuardTiles l_mode l_rank l_color_vectors l_reg_GTs
        l_guards = map (pShowAutoGuardString " && ") l_guardTiles
        l_guardNames = map (pSubstitute "!" "_Not_" . gName . fst) l_guardTiles
        l_tiles = map (pShowAutoTileString l_mode l_stencil) l_guardTiles
        l_str_GTs = zipWith (++) l_guards l_tiles
    in  (breakline ++ pShowHeader ++ 
            breakline ++ pShowColorVectors l_color_vectors ++ 
            breakline ++ concat l_str_GTs, 
            l_guardTiles)

pShowColorVectors :: [Homogeneity] -> String 
pShowColorVectors l_color_vectors = 
    breakline ++ "/* " ++ 
    breakline ++ show l_color_vectors ++ 
    breakline ++ " */" ++ breakline

pGetGuardTiles :: PMode -> Int -> [Homogeneity] -> [(PGuard, PTile)] -> [(PGuard, [PTile])]
pGetGuardTiles _ _ [] _ = []
pGetGuardTiles l_mode l_rank cL@(c:cs) l_reg_GTs =
    let l_len = size c
        l_gts = pGetGuardTilesTerm l_mode l_rank c l_reg_GTs  
    in  [l_gts] ++ pGetGuardTiles l_mode l_rank cs l_reg_GTs

pGetGuardTilesTerm :: PMode -> Int -> Homogeneity -> [(PGuard, PTile)] -> (PGuard, [PTile])
pGetGuardTilesTerm l_mode l_rank l_color l_reg_GTs =
    let l_len = size l_color
        l_guards = map fst l_reg_GTs
        l_guards' = zipWith pFillGuardOrder [0..] l_guards
        l_out_guard = foldr (pColorGuard l_rank l_color) emptyGuard l_guards'
        l_tiles = map snd l_reg_GTs
        l_tiles' = zipWith pFillTileOrder [0..] l_tiles
        l_out_tiles = foldr (pColorTile l_rank l_color) [emptyTile] l_tiles'
    in  (l_out_guard, l_out_tiles)

pColorGuard :: Int -> Homogeneity -> PGuard -> PGuard -> PGuard
pColorGuard l_rank l_color l_old_guard_l l_old_guard_r =
    let l_idx = gOrder l_old_guard_l
        l_condStr = 
            if testBit (o l_color) l_idx && testBit (a l_color) l_idx
               then [gName l_old_guard_l] ++ gComment l_old_guard_r 
               else if (not $ testBit (o l_color) l_idx) && 
                        (not $ testBit (a l_color) l_idx)
                       then ["!" ++ gName l_old_guard_l] ++ gComment l_old_guard_r
                       else gComment l_old_guard_r
    in  PGuard { gName = pSys $ intercalate "_" l_condStr, gRank = l_rank, gFunc = emptyGuardFunc, gComment = l_condStr, gOrder = 0 }
                            
pColorTile :: Int -> Homogeneity -> PTile -> [PTile] -> [PTile]
pColorTile l_rank l_color l_old_tile_l l_old_tiles_r =
    let l_idx = tOrder l_old_tile_l
        l_tiles = 
            if testBit (o l_color) l_idx && testBit (a l_color) l_idx
               then let l_new_tile_l = pSetTileOp PSERIAL l_old_tile_l
                    in  [l_new_tile_l] ++ l_old_tiles_r
               else if (testBit (o l_color) l_idx) &&
                        (not $ testBit (a l_color) l_idx)
                       then let l_new_tile_l = pSetTileOp PINCLUSIVE l_old_tile_l
                            in  [l_new_tile_l] ++ l_old_tiles_r
                       else l_old_tiles_r
    in  l_tiles

pShowHeader :: String
pShowHeader = 
    breakline ++ "#include <cstdio>" ++
    breakline ++ "#include <cstdlib>" ++
    breakline ++ "#include <cassert>" ++
    breakline ++ "#include <functional>" ++
    breakline ++ "#include <dlfcn.h>" ++ breakline 
