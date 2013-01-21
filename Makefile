PP_FILE=PMain.hs PMainParser.hs PBasicParser.hs PData.hs PShow.hs PUtils.hs PParser2.hs
pochoir : ${PP_FILE} 
	ghc -o pochoir -O --make PMain.hs
clean: 
	rm *.o *.hi 
