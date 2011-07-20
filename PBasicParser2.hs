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
    -- convert "Register_Kernel(g, k, ... ks) " to "Register_Obase_Kernel(g, k, bk)"
        do try $ pMember "Register_Kernel"
           (l_guard, l_kernels) <- parens pStencilRegisterKernelParams
           semi
           case Map.lookup l_id $ pStencil l_state of
               Nothing -> return (l_id ++ ".Register_Kernel(" ++ l_guard ++ ", " ++ 
                                  intercalate ", " l_kernels ++ 
                                  "); /* Register_Kernel with UNKNOWN Stencil " ++ 
                                  l_id ++ "*/" ++ breakline)
               Just l_stencil ->
                   do let l_pKernels = map (getValidKernel l_state) l_kernels
                      let l_pGuard = getValidGuard l_state l_guard
                      let l_validKernel = foldr (&&) True $ map fst l_pKernels
                      let l_validGuard = fst l_pGuard
                      if (l_validKernel == False || l_validGuard == False) 
                         then return (l_id ++ ".Register_Kernel(" ++ 
                                      l_guard ++ ", " ++ intercalate ", " l_kernels ++ 
                                      ");" ++ "/* Not all kernels are valid */ " ++ 
                                      breakline)
                         else do let l_regKernels = pShowRegKernel (pMode l_state) l_stencil (snd l_pGuard, map snd l_pKernels)
                                 return (l_regKernels)
    <|> do return (l_id)

