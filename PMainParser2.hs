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

-- The Main Parser for a second pass, the purpose of which is to transform
-- the user's specification based on the info collected in the first pass
module PMainParser2 where

import Text.ParserCombinators.Parsec

import Control.Monad

import PBasicParser
import PBasicParser2
import PUtils
import PData
import PShow
-- import Text.Show
import Data.List
import qualified Data.Map as Map

pToken1 :: GenParser Char ParserState String
pToken1 = 
        try pParsePochoirStencilMember1
--    <|> try pParsePochoirKernel1
    <|> do ch <- anyChar
           return [ch]
    <?> "line"

pParsePochoirStencilMember1 :: GenParser Char ParserState String
pParsePochoirStencilMember1 =
    do l_id <- try (pIdentifier)
       l_state <- getState
       try $ ppStencil1 l_id l_state

