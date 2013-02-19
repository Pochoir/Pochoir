pShowAutoTileString :: PMode -> PStencil -> (PName, [PTile]) -> String
pShowAutoTileString l_mode l_stencil (l_guardName, l_tiles@(t:ts)) =
    let l_comments = pShowAutoTileComments l_tiles
        -- getTileKernels also fills the guardFunc/tile_op into PKernelFuncs
        l_kernels = concatMap getTileKernels l_tiles
        l_unroll = foldr max 0 $ map pTileLength l_tiles
        -- group kernels by the tile index
        l_kernels_by_tIndex = groupBy eqIndexPKernel l_kernels
        l_tile_indices = map (kIndex . head) l_kernels_by_tIndex 
        -- l_kernel_funcs : [[PKernelFunc]]
        l_kernel_funcs = map (map kFunc) l_kernels_by_tIndex
        l_stmts = concatMap kfStmt $ concat l_kernel_funcs
        l_params = foldr union (kfParams $ head $ head l_kernel_funcs) 
                               (map kfParams $ tail $ head l_kernel_funcs)
        l_iters = getIterFromKernel l_mode l_stencil l_params l_stmts
        -- l_rev_kernel_funcs : [[PKernelFunc]]
        l_rev_kernel_funcs = map (map (pFillIters l_iters)) l_kernel_funcs
        l_id = sName l_stencil
    in  case l_mode of
            PAllCondTileMacroOverlap ->
                pSplitAllCondTileOverlapScope
                  ("Macro_", l_mode, l_id, l_guardName, l_unroll, l_tile_indices, l_rev_kernel_funcs, l_stencil, l_comments)
                  (pShowAllCondTileOverlapMacroKernels False)

pSplitAllCondTileOverlapScope :: (String, PMode, PName, PName, Int, [[Int]], [[PKernelFunc]], PStencil, String) -> (String -> PStencil -> [[Int]] -> [[PKernelFunc]] -> String) -> String
pSplitAllCondTileOverlapScope (l_tag, l_mode, l_id, l_guardName, l_unroll, l_tile_indices, l_kfss, l_stencil, l_comments) l_showAllCondTileOverlapKernel = 
    let oldKernelName = intercalate "_" $ map kfName $ concat l_kfss
        bdryKernelName = l_tag ++ "boundary_" ++ oldKernelName
        obaseKernelName = l_tag ++ "interior_" ++ oldKernelName
        regBound = sRegBound l_stencil
        l_pShape = pSysShape $ foldr mergePShapes emptyShape (map kfShape $ concat l_kfss)
        bdryKernel = if regBound
                        then pShowAllCondTileOverlapMacroKernels True bdryKernelName
                                l_stencil l_pShape l_tile_indices l_kfss
                        else ""
        -- usually, the l_showAllCondTileOverlapKernel is just
        -- (pShowAllCondTileOverlapMacroKernels False)
        obaseKernel = l_showAllCondTileOverlapKernel obaseKernelName
                        l_stencil l_pShape l_tile_indices l_kfss
        runKernel = if regBound
                       then obaseKernelName ++ ", " ++ bdryKernelName
                       else obaseKernelName
    in  (breakline ++ l_comments ++
         breakline ++ show l_pShape ++ bdryKernel ++ breakline ++ obaseKernel ++
         breakline ++ l_id ++ ".Register_Tile_Obase_Kernels(" ++ l_guardName ++
         ", " ++ show l_unroll ++ ", " ++ runKernel ++ ");" ++ breakline)

pShowAllCondTileOverlapMacroKernels :: Bool -> String -> PStencil -> PShape -> [[Int]] -> [[PKernelFunc]] -> String
pShowAllCondTileOverlapMacroKernels l_bound l_name l_stencil l_pShape l_tile_indices l_kfss@(k:ks) =
    let l_t = "t"
        l_t_begin = l_t ++ "0"
        l_t_end = l_t ++ "1"
        l_rank = sRank l_stencil
        l_arrayInUse = sArrayInUse l_stencil
        -- We are assuming that all kernel functions have the 
        -- same number of input parameters
        l_kfParams = (kfParams . head . head) l_kfss
        l_defMacro = if l_bound
                        then pDefMacroArrayInUse "boundary" l_arrayInUse l_kfParams ++
                                breakline ++ pDefPMODLU
                        else pDefMacroArrayInUse "interior" l_arrayInUse l_kfParams
        l_undefMacro = if l_bound
                          then pUndefPMODLU ++ breakline ++
                                pUndefMacroArrayInUse l_arrayInUse l_kfParams
                          else pUndefMacroArrayInUse l_arrayInUse l_kfParams
        l_showPhysGrid = "Grid_Info <" ++ show l_rank ++ "> l_phys_grid = " ++
                            sName l_stencil ++ ".get_phys_grid();"
        l_unfold_kernels = pShowAllCondTileOverlapSingleMacroKernel l_t l_tile_indices l_kfss
        l_kernelFuncName = pSys l_name
        l_header = "/* KNOWN! */ auto " ++ l_kernelFuncName ++
                   " = [&] (int " ++ l_t_begin ++ ", int " ++ l_t_end ++ ", " ++ 
                   " Grid_Info <" ++ show l_rank ++ "> const & grid) {"
        l_tail = "};" ++ breakline ++ "Pochoir_Obase_Kernel <" ++ show l_rank ++
                 "> " ++ l_name ++ "( " ++ l_kernelFuncName ++ ", " ++
                 shapeName l_pShape ++ " );"
    in  breakline ++ l_defMacro ++
        breakline ++ l_header ++
        breakline ++ "Grid_Info <" ++ show l_rank ++ "> l_grid = grid;" ++
        breakline ++ l_showPhysGrid ++
        breakline ++ pShowTimeLoopHeader l_t l_t_begin l_t_end ++
        breakline ++ pShowMetaGridHeader l_bound (tail l_kfParams) ++
        l_unfold_kernels ++
        breakline ++ pShowMetaGridTail (tail l_kfParams) ++
        breakline ++ pAdjustTrape l_rank ++
        breakline ++ pShowTimeLoopTail ++
        breakline ++ l_tail ++
        breakline ++ l_undefMacro ++ breakline

pShowAllCondTileOverlapSingleMacroKernel :: String -> [[Int]] -> [[PKernelFunc]] -> String
pShowAllCondTileOverlapSingleMacroKernel _ [] _ = ""
pShowAllCondTileOverlapSingleMacroKernel l_t l_tile_indices@(t:ts) l_kL@(k:ks) = 
    let l_spatial_params = tail $ kfParams $ head k
        l_rank = length l_spatial_params
        l_params = [l_t] ++ l_spatial_params
        l_dim_sizes = getTileSizes l_tile_indices
        l_tile_index = head l_tile_indices
        l_guard_head = pShowTileGuardHeadOnAll l_params l_dim_sizes l_tile_index
                                                (shapeTimeShift $ kfShape $ head k)
        l_guard_tail = pShowUnrollGuardTail $ length ts
        k' = groupBy eqTileOpPKernelFunc k
        l_kernels = pShowOverlapMacroKernels l_t l_spatial_params k'
    in  breakline ++ l_guard_head ++
        breakline ++ l_kernels ++
        breakline ++ l_guard_tail ++
        pShowAllCondTileOverlapSingleMacroKernel l_t ts ks

-- pShowOverlapMacroKernels will show the macro kernels with the same tile_index
-- All kernel functions are groupBy kfTileOp
pShowOverlapMacroKernels :: String -> [PName] -> [[PKernelFunc]] -> String
pShowOverlapMacroKernels _ _ [] = ""
pShowOverlapMacroKernels l_t l_spatial_params l_kfss@(k:ks) =
    case (kfTileOp . head) k of
        PSerial -> concatMap (show . kfStmt) k ++ 
                   pShowOverlapMacroKernels l_t l_spatial_params ks
        PEXCLUSIVE -> pShowInclusiveMacroKernels PEXCLUSIVE l_t l_spatial_params k ++
                      pShowOverlapMacroKernels l_t l_spatial_params ks
        PINCLUSIVE ->  pShowInclusiveMacroKernels PINCLUSIVE l_t l_spatial_params k ++
                      pShowOverlapMacroKernels l_t l_spatial_params ks

pShowInclusiveMacroKernels :: TileOp -> String -> [PName] -> [PKernelFunc] -> String
pShowInclusiveMacroKernels _ _ _ [] = ""
pShowInclusiveMacroKernels l_tile_op l_t l_spatial_params l_kfs@(k:ks) =
    let g = (gName $ kfGuardFunc k) ++ " ( " ++
            l_t ++ ", " ++ intercalate ", " l_spatial_params ++
            " ) "
        l_tail = if (l_tile_op == PEXCLUSIVE && length ks > 0) 
                    then "} else "
                    else "}"
    in  (breakline ++ "if (" ++ g ++ ") {" ++ 
         breakline ++ (show $ kfStmt k) ++ 
         breakline ++ l_tail ++ breakline ++ 
         pShowOverlapMacroKernels l_t l_spatial_params ks)

