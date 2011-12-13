PP_FILE=PMain.hs PMainParser.hs PBasicParser.hs PData.hs PShow.hs PUtils.hs PMainParser2.hs PBasicParser2.hs
GEN_FILE=PGenMain.hs PGenMainParser.hs PGenBasicParser.hs PGenData.hs PGenShow.hs PGenUtils.hs PGenMainParser2.hs PGenBasicParser2.hs
pochoir : ${PP_FILE} 
	ghc -o pochoir -O --make PMain.hs
genkernels : ${GEN_FILE} 
	ghc -o genkernels -O --make PGenMain.hs
clean: 
	rm *.o *.hi pochoir 
