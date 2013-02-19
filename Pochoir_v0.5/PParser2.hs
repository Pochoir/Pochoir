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

module PParser2 where
{- modules for 2nd pass parsers -}

import Text.ParserCombinators.Parsec
import Control.Monad

import PBasicParser
import PShow
import PData
import qualified Data.Map as Map

pToken1 :: GenParser Char ParserState String
pToken1 = do reserved "Pochoir_Array"
             (l_type, l_rank) <- angles pDeclStatic <?> "Pochoir_Array static parameters"
             l_arrayDecl <- commaSep1 pDeclDynamic <?> "Pochoir_Array Dynamic parameters"
             semi 
             l_state <- getState
             return (concat $ 
                   map (outputSArray (l_type, l_rank) (pArray l_state)) l_arrayDecl)
             -- return ("Pochoir_Array <" ++ show l_type ++ ", " ++ show l_rank ++ "> " ++ pShowArrayDynamicDecl l_arrayDecl ++ ";\n")
      <|> do ch <- anyChar
             return [ch]
      <?> "line"

outputSArray :: (PType, Int) -> Map.Map PName PArray -> ([PName], PName, [DimExpr]) -> String
outputSArray (l_type, l_rank) m_array (l_qualifiers, l_array, l_dims) =
    case Map.lookup l_array m_array of
        Nothing -> breakline ++ "Pochoir_Array <" ++ 
                    show l_type ++ ", " ++ show l_rank ++ "> " ++ 
                    pShowDynamicDecl [(l_qualifiers, l_array, l_dims)] pShowArrayDim ++ ";" 
        Just l_pArray -> breakline ++ "Pochoir_Array <" ++ 
                    show l_type ++ ", " ++ show l_rank ++ "> " ++ 
                    pShowDynamicDecl [(l_qualifiers, l_array, l_dims)] pShowArrayDim ++ ";" 


