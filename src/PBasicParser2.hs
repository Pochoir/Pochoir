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

-- The Basic Parser for a second pass, the purpose of which is to transform
-- the user's specification based on the info collected in the first pass
module PBasicParser2 where

import Text.ParserCombinators.Parsec
import qualified Text.ParserCombinators.Parsec.Token as Token
import Text.ParserCombinators.Parsec.Expr
import Text.ParserCombinators.Parsec.Language
import Text.Read (read)

import Data.Char
import Data.List
import qualified Data.Map as Map

import PBasicParser
import PShow
import PUtils
import PData

ppStencil1 :: String -> ParserState -> GenParser Char ParserState String
ppStencil1 l_id l_state = 
        do try $ pMember "Run"
           l_tstep <- parens exprStmtDim
           semi
           let l_mode = pMode l_state
           case Map.lookup l_id $ pStencil l_state of
               Nothing -> return (l_id ++ ".Run(" ++ show l_tstep ++ ");")
               Just l_stencil -> 
                   do let l_inFile = pInFile l_state
                      let l_fname = pSubstitute ".cpp" "" l_inFile
                      let l_stencil_name = sName l_stencil
                      let l_arrayInUse = map aName $ sArrayInUse l_stencil
                      let l_runCmd = if l_mode == PMUnroll 
                            || l_mode == PCondMacro
                            || l_mode == PCondCPointer
                            || l_mode == PCondPointer
                            || l_mode == PCondOptPointer
                          then ".Run_Obase_All_Cond(" 
                          else ".Run_Obase_Unroll_T(" 
                      return (l_id ++ l_runCmd ++ show l_tstep ++
                              ", " ++ (mkQuote . show) l_mode ++
                              ", " ++ mkQuote l_fname ++
                              ", " ++ mkQuote l_stencil_name ++
                              ", " ++ (intercalate ", " l_arrayInUse) ++
                              "); /* KNOWN */" ++ breakline)
    <|> do return (l_id)

