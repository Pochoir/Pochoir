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

module PGenData where

import Text.ParserCombinators.Parsec
import Control.Monad

import Data.Char
import Data.List
import Data.Bits
import qualified Data.Map as Map

type PName = String
type PValue = Int
type Uop = String 
type Bop = String
breakline :: String
breakline = "\n\t"

stripWhite :: String -> String
stripWhite l = dropWhile isSpace l

externC :: String
externC = "extern \"C\" " ++ breakline

pTab :: String
pTab = "\t"

-- We use newtype just because we want to manually derive a Show instance of PShift
newtype PShift = PShift Int deriving (Eq, Ord)
data RegionT = Periodic | Nonperiodic | UnknownRegionT deriving Show
data PBasicType = PChar | PShort | PLong | PSigned | PUnsigned | PVoid | PInt | PDouble | PFloat | PBool | PUserType deriving Eq
data PPred = PLEQ | PGEQ | PEQ | PNEQ deriving (Eq, Show)
data PRWMode = PRead | PWrite | POther deriving (Eq, Show)
data PType = PType {
    basicType :: PBasicType,
    typeName :: String
} deriving Eq
data PMode = PHelp | PDefault | PDebug | PCaching | PCPointer | PMUnroll | POptPointer | PPointer | PMacroShadow | PNoPP | PAllCondTileMacro | PAllCondTileCPointer | PAllCondTilePointer | PAllCondTileOptPointer | PUnrollTimeTileMacro | PUnrollTimeTileCPointer | PUnrollTimeTilePointer | PUnrollTimeTileOptPointer | PAllCondTileMacroOverlap | PAllCondTileCPointerOverlap | PAllCondTilePointerOverlap | PAllCondTileOptPointerOverlap | PUnrollTimeTileMacroOverlap | PUnrollTimeTileCPointerOverlap | PUnrollTimeTilePointerOverlap | PUnrollTimeTileOptPointerOverlap deriving Eq
data TileOp = PNULL | PSERIAL | PEXCLUSIVE | PINCLUSIVE deriving (Eq, Show)

data Homogeneity = Homogeneity {
    size :: Int, -- width of the color vector
    o :: Int, -- bit-wise 'or' of the color vector
    a :: Int  -- bit-wise 'and' of the color vector
} deriving Eq

data PArray = PArray {
    aName :: PName,
    aType :: PType,
    aRank :: Int,
    aMaxShift :: Int,
    aToggle :: Int,
    aDims :: [DimExpr],
    aRegBound :: Bool,
    aComment :: String
} deriving (Show, Eq)

data PStencil = PStencil {
    sName :: PName,
    sRank :: Int,
    sToggle :: Int,
    sTimeShift :: Int,
    sArrayInUse :: [PArray],
    sShape :: PShape,
    sRegTileKernel :: [(PGuard, PTile)], -- the original user registered (guard, tile) pair
    sRegBound :: Bool,
    sComment :: String
} deriving Show

data PShape = PShape {
    shapeName :: PName,
    shapeRank :: Int,
    shapeLen :: Int,
    shapeToggle :: Int,
    shapeSlopes :: [Int],
    shapeTimeShift :: Int,
    shape :: [[Int]],
    shapeComment :: String
} deriving Eq

data PRange = PRange {
    rName :: PName,
    rFirst :: DimExpr,
    rLast :: DimExpr,
    rStride :: DimExpr,
    rComment :: String
} deriving Show 

-- (NameOfIter, arrayInUse, correspondingDimExpr, RWMode)
type Iter = (String, PArray, [DimExpr], PRWMode)

data PKernelFunc = PKernelFunc {
    kfName :: PName,
    kfParams :: [PName],
    kfStmt :: [Stmt],
    -- How many lines of statement do we have
    kfStmtSize :: Int,
    kfIter :: [Iter],
    kfShape :: PShape,
    kfTileOp :: TileOp,
    kfTileOrder :: Int,
    kfTileSizes :: [Int], -- sizes of each time/space dimension in original tile
    kfTileIndex :: [Int], -- index in the tile
    kfGuardFunc :: PGuardFunc, -- each kernel function is guarded by one and only one 
                               -- guard function
    kfComment :: String
} deriving Show

data PGuardFunc = PGuardFunc {
    gfName :: PName,
    gfParams :: [PName],
    gfStmt :: [Stmt],
    -- How many lines of statement do we have
    gfStmtSize :: Int,
    gfIter :: [Iter],
    gfOrder :: Int,
    gfComment :: String
} deriving Show

data PGuard = PGuard {
    gName :: PName,
    gRank :: Int,
    gFunc :: PGuardFunc,
    gColor :: Homogeneity,
    gOrder :: Int,
    gComment :: [String]
} deriving Show

data PKernel = PKernel {
    kName :: PName,
    kRank :: Int,
    kFunc :: PKernelFunc,
    kShape :: PShape,
    kIndex :: [Int], -- this is the index in the tile
    kTileSizes :: [Int], -- The sizes of each time/space dimension in original tile
    kComment :: String
} deriving Show

data PTileKernel = SK PKernel 
                 | LK [PTileKernel]

data PTile = PTile {
    tName :: PName,
    tRank :: Int,
    tSizes :: [Int], -- size of each time/spatial dimension
    -- tKernel :: PTileKernel,
    tKernels :: [PKernel],
    tOrigGuard :: PGuard, -- each tile has one and only one guard
    tOp :: TileOp, -- the OP is used in determining how to tile in automatically generated overlapped kernel
    tOrder :: Int, -- tile order
    tColor :: Homogeneity,
    tComment :: String
} deriving Show

data PMTile = PMTile {
    -- mtSizes == concatMap mttSizes mtItems
    mtSizes :: [Int], -- sizes of each space/time dimensions, starting from t, i, j, ... 
    mtLatestTileOrder :: Int, -- the tile-order of latest inserted kernel_func
    mtKernelFuncs :: [PKernelFunc], -- collection of all inserted kernel funcs
    mtTerms :: [PMTileTerm]
} 

data PMTileTerm = PMTileTerm {
    mttIndex :: [Int],
    mttSizes :: [Int],
    mttLatestTileOrder :: Int, -- the tile-order of latest inserted kernel_func
    mttItem :: PMTileItem 
} 
    
data PMTileItem = ST [PKernelFunc]
                | MT [PMTileTerm]

emptyPType :: PType
emptyPType = PType { basicType = PUserType, typeName = "" }

emptyPArray :: PArray
emptyPArray = PArray { aName = "", aType = emptyPType, aRank = 0, aMaxShift = 0, aToggle = 0, aDims = [], aRegBound = False, aComment = cEmpty "Pochoir_Array" }

emptyShape :: PShape
emptyShape = PShape { shapeName = "", shapeRank = 0, shapeLen = 0, shapeToggle = 0, shapeSlopes = [], shapeTimeShift = 0, shape = [], shapeComment = cEmpty "Pochoir_Shape" }

emptyKernelFunc :: PKernelFunc
emptyKernelFunc = PKernelFunc { kfName = "", kfParams = [], kfStmt = [], kfStmtSize = 0, kfIter = [], kfShape = emptyShape, kfTileOp = PNULL, kfGuardFunc = emptyGuardFunc, kfTileOrder = 0,kfTileSizes = [], kfTileIndex = [], kfComment = cEmpty "Pochoir_Kernel_Func" }

emptyKernel :: PKernel
emptyKernel = PKernel { kName = "", kRank = 0, kFunc = emptyKernelFunc, kShape = emptyShape, kIndex = [], kTileSizes = [], kComment = cEmpty "Pochoir_Kernel" }

emptyGuardFunc :: PGuardFunc
emptyGuardFunc = PGuardFunc { gfName = "", gfParams = [], gfStmt = [], gfStmtSize = 0, gfIter = [], gfOrder = 0, gfComment = cEmpty "Pochoir_Guard_Func" }

emptyGuard :: PGuard
emptyGuard = PGuard { gName = "", gRank = 0, gFunc = emptyGuardFunc, gComment = [], gOrder = 0, gColor = emptyColor }

emptyColor :: Homogeneity
emptyColor = Homogeneity { size = 0, o = 0, a = 0 }

emptyTileKernel :: PTileKernel
emptyTileKernel = LK [] 

emptyTile :: PTile
emptyTile = PTile { tName = "", tRank = 0, tSizes = [], tKernels = [], tComment = cEmpty "Pochoir_Tile", tOp = PNULL, tOrigGuard = emptyGuard, tOrder = 0, tColor = emptyColor }

-- prefix 'c' means "comment"
cUnknown :: String -> String
cUnknown l_name = "/* Unknown " ++ l_name ++ " */"

-- prefix 'c' means "comment"
cKnown :: String -> String
cKnown l_name = "/* Known " ++ l_name ++ " */"

-- prefix 'c' means "comment"
cEmpty :: String -> String
cEmpty l_name = "/* Empty " ++ l_name ++ " */"

data ParserState = ParserState {
    pMode  :: PMode,
    pColorVectors :: [Homogeneity],
    pArray :: Map.Map PName PArray,
    pStencil :: Map.Map PName PStencil,
    pRange :: Map.Map PName PRange,
    pShape :: Map.Map PName PShape,
    pKernelFunc :: Map.Map PName PKernelFunc,
    pKernel :: Map.Map PName PKernel,
    pGuardFunc :: Map.Map PName PGuardFunc,
    pGuard :: Map.Map PName PGuard,
    pTile :: Map.Map PName PTile,
    pTileOrder :: Int,
    pGuardOrder :: Int,
    pGuardFuncOrder :: Int,
    pGenPlan :: Map.Map Int PStencil, -- this is a snapshot of each PStencil object
                                      -- at the time of calling Gen_Plan_Obase(T)
    pGenPlanOrder :: Int,
    pColorNum :: Int
} deriving Show

data Expr = VAR String String 
          -- PVAR : (*&)p(a, b, ...)
          | PVAR String String [DimExpr]
          -- BVAR : v[a]
          | BVAR String DimExpr
          -- BExprVAR : v[a+b+c]
          | BExprVAR String Expr
          -- SVAR : type(p(a, b)).field
          | SVAR PType Expr String String
          -- PSVAR : (type*)(p(a, b)).field
          | PSVAR PType Expr String String
          -- Uno is prefix unary operator
          | Uno Uop Expr
          -- PostUno is suffix unary operator
          | PostUno Uop Expr
          | Duo Bop Expr Expr
          | INT Int
          | FLOAT Double
          | BOOL String
          | PARENS Expr
          deriving Eq

data Stmt = BRACES [Stmt]
          | EXPR Expr
          | DEXPR [PName] PType [Expr]
          | IF Expr Stmt Stmt
          | SWITCH Expr [Stmt]
          | CASE PValue [Stmt]
          | DEFAULT [Stmt]
          | NOP
          | BREAK
          | DO Expr [Stmt]
          | WHILE Expr [Stmt]
          | FOR [[Stmt]] Stmt
          | CONT
          | RET Expr
          | RETURN
          | UNKNOWN String
          deriving Eq

data DimExpr = DimVAR String 
          | DimDuo Bop DimExpr DimExpr
          | DimParen DimExpr
          | DimINT Int
          deriving Eq

instance Ord DimExpr where
    (<) a b = (<) (transDimExprFloat a) (transDimExprFloat b)
    (>) a b = (>) (transDimExprFloat a) (transDimExprFloat b)
    (<=) a b = (<=) (transDimExprFloat a) (transDimExprFloat b)
    (>=) a b = (>=) (transDimExprFloat a) (transDimExprFloat b)
    compare a b = compare (transDimExprFloat a) (transDimExprFloat b)

transDimExprFloat :: DimExpr -> Float
transDimExprFloat (DimVAR v) = 1
transDimExprFloat (DimDuo "+" a b) = (transDimExprFloat a) + (transDimExprFloat b)
transDimExprFloat (DimDuo "-" a b) = (transDimExprFloat a) - (transDimExprFloat b)
transDimExprFloat (DimDuo "*" a b) = (transDimExprFloat a) * (transDimExprFloat b)
transDimExprFloat (DimDuo "/" a b) = (transDimExprFloat a) / (transDimExprFloat b)
transDimExprFloat (DimParen e) = (transDimExprFloat e)
transDimExprFloat (DimINT n) = fromIntegral n
transDimExprFloat _ = 0

instance Show PMode where
    show PHelp = "-help" 
    show PDefault = "-default" 
    show PDebug = "-debug" 
    show PCaching = "-split-caching" 
    show PCPointer = "-split-c-pointer" 
    show POptPointer = "-split-opt-pointer" 
    show PPointer = "-split-pointer" 
    show PMacroShadow = "-split-macro-shadow" 
    show PNoPP = "-No-Preprocessing"
    show PMUnroll = "-unroll-multi-kernel"
    show PAllCondTileMacro = "-all-cond-tile-macro"
    show PAllCondTileCPointer = "-all-cond-tile-c-pointer"
    show PAllCondTilePointer = "-all-cond-tile-pointer"
    show PAllCondTileOptPointer = "-all-cond-tile-opt-pointer"
    show PUnrollTimeTileMacro = "-unroll-t-tile-macro"
    show PUnrollTimeTileCPointer = "-unroll-t-tile-c-pointer"
    show PUnrollTimeTilePointer = "-unroll-t-tile-pointer"
    show PUnrollTimeTileOptPointer = "-unroll-t-tile-opt-pointer"
    show PAllCondTileMacroOverlap = "-all-cond-tile-macro-overlap"
    show PAllCondTileCPointerOverlap = "-all-cond-tile-c-pointer-overlap"
    show PAllCondTilePointerOverlap = "-all-cond-tile-pointer-overlap"
    show PAllCondTileOptPointerOverlap = "-all-cond-tile-opt-pointer-overlap"
    show PUnrollTimeTileMacroOverlap = "-unroll-t-tile-macro-overlap"
    show PUnrollTimeTileCPointerOverlap = "-unroll-t-tile-c-pointer-overlap"
    show PUnrollTimeTilePointerOverlap = "-unroll-t-tile-pointer-overlap"
    show PUnrollTimeTileOptPointerOverlap = "-unroll-t-tile-opt-pointer-overlap"

instance Show Homogeneity where
    show h = "<" ++ showBin (size h - 1) (o h) ++ ", " ++ showBin (size h - 1) (a h) ++ ">\n"

{-
instance Show PMTile where
    show a = let l_sizes = mtSizes a
                 l_terms = mtTerms a
             in  breakline ++ "PMTile (" ++ 
                 "sizes = " ++ show l_sizes ++ 
                 "), terms = " ++ breakline ++ show l_terms

instance Show PMTileTerm where
    show a = let l_indices = mttIndex a
                 l_sizes = mttSizes a
                 l_item = mttItem a
             in  pTab ++ "PMTileTerm (" ++ 
                 "indices = " ++ show l_indices ++ 
                 ", sizes = " ++ show l_sizes ++ 
                 "), item = " ++ show l_item

instance Show PMTileItem where
    show (ST l_ks) = concatMap kfName l_ks ++ breakline
    show (MT l_mts) = breakline ++ pTab ++ show l_mts
 -}

pShowPMTile :: Int -> PMTile -> String
pShowPMTile l_rec l_mtile =
    let l_sizes = mtSizes l_mtile
        l_terms = mtTerms l_mtile
        l_tile_order = mtLatestTileOrder l_mtile
    in  breakline ++ "PMTile (sizes = " ++
        show l_sizes ++ ", latest_tile_order = " ++ show l_tile_order ++
        "), terms = " ++ breakline ++ 
        (concatMap (pShowPMTileTerm (l_rec + 1)) l_terms)

pShowPMTileTerm :: Int -> PMTileTerm -> String
pShowPMTileTerm l_rec l_mterm =
    let l_indices = mttIndex l_mterm
        l_sizes = mttSizes l_mterm
        l_tile_order = mttLatestTileOrder l_mterm
        l_item = mttItem l_mterm
    in  (concat $ replicate l_rec pTab) ++
        "PMTileTerm (indices = " ++ show l_indices ++
        ", sizes = " ++ show l_sizes ++
        ", latest_tile_order = " ++ show l_tile_order ++
        "), item = " ++ pShowPMTileItem (l_rec + 1) l_item

pShowPMTileItem :: Int -> PMTileItem -> String
pShowPMTileItem l_rec (ST l_ks) = (intercalate ", " $ map kfName l_ks) ++ breakline
pShowPMTileItem l_rec (MT l_mts) = 
    breakline ++ 
    (concatMap (pShowPMTileTerm l_rec) l_mts)

showBin :: Int -> Int -> String
showBin (-1) a = ""
showBin n a = showBinTerm n a ++ showBin (n-1) a
    where showBinTerm n a = if testBit a n then "1" else "0"

instance Show PType where
    show ptype = typeName ptype

instance Show PShape where
    show l_pShape = let l_rank = shapeRank l_pShape
                        l_name = shapeName l_pShape
                        l_shapes = shape l_pShape
                        l_toggle = shapeToggle l_pShape
                        l_slopes = shapeSlopes l_pShape
                    in  breakline ++ "Pochoir_Shape <" ++ show l_rank ++ 
                        "> " ++ l_name ++ " [ ] = " ++ pShowShapes l_shapes ++ 
                        ";" ++ breakline ++ "/* toggle: " ++ show l_toggle ++ 
                        "; slopes: " ++ show l_slopes ++ " */" ++ breakline

instance Show PTileKernel where
    show (SK k) = kName k
    show (LK kList) = "{" ++ showList kList "" ++ "}"
    showList [] = showString ""
    showList (k:ks) = shows k . showL ks
        where showL [] = showString ""
              showL (k:ks) = showString ", " . shows k . showL ks

instance Show Expr where
    show (VAR q str) = q ++ str
    show (PVAR q a xList@(x:xs)) = q ++ a ++ "(" ++ showList xList "" ++ ")"
    show (BVAR a x) = a ++ "[" ++ show x ++ "]"
    show (BExprVAR a e) = a ++ "[" ++ show e ++ "]"
    show (SVAR t e c f) = show t ++ "(" ++ show e ++ ")" ++ c ++ f
    show (PSVAR t e c f) = "(" ++ show t ++ ")(" ++ show e ++ ")" ++ c ++ f
    show (Uno uop expr) = uop ++ show expr 
    show (PostUno uop expr) = show expr ++ uop
    show (Duo bop lexpr rexpr) = show lexpr ++ " " ++ bop ++ " " ++ show rexpr
    show (PARENS expr) = "(" ++ show expr ++ ")"
    show (INT n) = show n
    show (FLOAT n) = show n
    show (BOOL b) = b
    showList [] = showString ""
    showList (x:xs) = showString breakline . shows x . showChar ';' . showList xs

instance Show DimExpr where
    show (DimVAR str) = str
    show (DimDuo bop lexpr rexpr) = show lexpr ++ " " ++ bop ++ " " ++ show rexpr
    show (DimParen e) = "(" ++ show e ++ ")"
    show (DimINT n) = show n
    showList (x:xs) = shows x . showl xs
                      where showl [] = showString "" 
                            showl (x:xs) = showString ", " . shows x . showl xs

instance Show Stmt where
    show (BRACES tL@(t:ts)) = "{" ++ showList tL "" ++ breakline ++ "}"
    show (EXPR expr) = show expr ++ ";"
    show (DEXPR qs declType es) = intercalate " " qs ++ " " ++ 
        (intercalate ", " $ map show es) ++ ";"
    show (IF expr l_stmt NOP) = "if " ++ show expr ++ 
        breakline ++ show l_stmt 
    show (IF expr l_stmt r_stmt) = "if (" ++ show expr ++ ") {" ++ 
        breakline ++ show l_stmt ++ 
        breakline ++ "} else {" ++ show r_stmt ++ "}"
    show (SWITCH expr tL@(t:ts)) = "switch " ++ show expr ++ "{" ++
        showList tL "" ++ breakline ++ "} /* end of switch */" ++ breakline
    show (CASE l_value tL@(t:ts)) = "case " ++ show l_value ++ " : " ++
        showList tL ""
    show (DEFAULT tL@(t:ts)) = "default :" ++ showList tL "" 
    show (NOP) = ""
    show (BREAK) = "break;"
    show (DO expr stmts) = "do {" ++ 
                            breakline ++ showList stmts "" ++ breakline ++ 
                            " }while " ++ show expr ++ ";" ++ breakline
    show (WHILE expr stmts) = "while (" ++ show expr ++ ") {" ++
        showList stmts "" ++ breakline ++ "} /* end of while */" ++ breakline
    show (FOR ttL@(t:ts) l_stmt) = "for " ++ showForListList ttL ++
        breakline ++ show l_stmt
            where showForListList (t:ts) = "(" ++ showForList t ++ showForListListL ts
                  showForListListL [] = ")"
                  showForListListL (t:ts) = "; " ++ showForList t ++ showForListListL ts
                  showForList (t:ts) = showForExpr t ++ showForListL ts
                  showForListL [] = ""
                  showForListL (t:ts) = ", " ++ showForExpr t ++ showForListL ts
                  showForExpr NOP = ""
                  showForExpr (EXPR expr) = show expr 
                  showForExpr (DEXPR qs declType expr) = intercalate " " qs ++ 
                                                         " " ++ show expr
    show (CONT) = "continue;" ++ breakline
    show (RET e) = "return (" ++ show e ++ ");" ++ breakline
    show (RETURN) = "return;" ++ breakline
    show (UNKNOWN stmt) = stmt ++ "\t/* Unrecognized! */" 
    showList [] = showString ""
    showList (x:xs) = shows x . showString breakline . showList xs
    -- showList (x:xs) = showString breakline . shows x . showList xs

instance Show PShift where
    show (PShift n) = show n
    showList [] = showString ""
    showList (x:xs) = showChar '{' . shows x . showl xs
                        where showl [] = showChar '}'
                              showl (x:xs) = showString ", " . shows x . showl xs

pShowShapes :: [[Int]] -> String
pShowShapes [] = ""
pShowShapes aL@(a:as) = "{" ++ pShowShape a ++ pShowShapesL as
    where pShowShapesL [] = "}"
          pShowShapesL (x:xs) = ", " ++ pShowShape x ++ pShowShapesL xs

pShowShape :: [Int] -> String
pShowShape [] = ""
pShowShape aL@(a:as) = "{" ++ show a ++ pShowShapeL as
    where pShowShapeL [] = "}"
          pShowShapeL (x:xs) = ", " ++ show x ++ pShowShapeL xs

