#Script to create color and shapes vectors for the AFC plots of Ecklonia SA
#Includes code for the vectors and population listing + color matching

#Turn Maxima pops into FALSE
maximaPops <- AFCdataInd$Population > 22

#Count the number of FALSES (nº of maxima samples)
length(maximaPops[maximaPops == FALSE])

#Counte the number of TRUES (nº of radiata samples)
length(maximaPops[maximaPops == TRUE])

#Create vector for 2 colors AFC, based on the FALSE and TRUE counts
k2colors <- rep(c("#0099E6", "#FF6600"), times = c(507, 266))

#Create vector for shapes according to sampling site/species
plotShapes <- rep(c(16,17,15,16,17,15), times = c(264,97,22,46,78,266))


#Create vector for 8 colors, to match the Structure plot. Counting was done by checking the column number of individuals
#and populations on the AFC table, plus using a randomvector to count the amount of samples in each group
k8colors <- rep(c("#2242D3","#0099E6", "#2EBFD2", "#0099E6", "#2EBFD2", "#6F65FA", "#15B79D", "#0099E6","#3CFFF6", "#EBB35C", "#FF6600"),
                times = c(96, 96, 24, 72, 25, 48, 23, 46, 77, 112, 154))




randomvector <- c(620:773)

##########################################################################################################################


Color per pop, according to hap network:
  pop = EMDP	1  #00AAFF
  pop = EMPN	2 #0000FF
  pop = EMKZ	3 #0000FF
  pop = EMHK	4 #0000FF
  pop = EMDB	5 #0000FF
  pop = EMJB	6 #0000FF #008FB7 (JB03)
  pop = EMKM	7 #0000FF
  pop = EMMP	8 #4CE5E5
  pop = EMVM	9 #89B236 #C6FF6A 07FF1B
  pop = EMSJ	10  #4CE5E5
  pop = EM-MUI	11 #A6B7FF #AAFFFF (01) #B482FF (16)
  pop = EMBB	12  #4CE5E5 #7389CE (09)
  pop = EM-BOT	13  #4CE5E5
  pop = ER-h-BOT	14 #FFD47A
  pop = EXBB	15 #FFD47A
  pop = ERRB	16 #00AA00
  pop = EMQP	17  #4CE5E5 A6B7FF
  pop = EMCA	18
  pop = EMDH	19 #A6B7FF
  pop = ER-HOO	20 #A6B7FF
  pop = ERDH	21 #A6B7FF #FFC800 (02)
  pop = ERSDH	22  #FF8C00
  pop = ERPE	23 #FF8C00
  pop = ERPA	24  #CF8A00
  pop = ERTS	25 #E56A0C
  pop = ERDM	26 #FF6600
  pop = ERDW	27 #FFC800
  pop = ERMK	28 #FFC800
  pop = ESTL	29 #F2E26C
  pop = ER-MOZ	30 #FFC800 
  pop = CBS	31#FFC800 FFF12A
  pop = E	32 #FFC800  FFF12A
  pop = PS	33 #FFC800  FFF12A
  pop = SH	34 #FFC800 FFF12A
   
  
  Color per pop, according to structure K=8:
    pop = EMDP	1  #2242D3
  pop = EMPN	2 #2242D3
  pop = EMKZ	3 #2242D3
  pop = EMHK	4 #2242D3
  pop = EMDB	5 #0099E6
  pop = EMJB	6 #0099E6
  pop = EMKM	7 #0099E6
  pop = EMMP	8 #0099E6
  pop = EMVM	9 #2EBFD2
  pop = EMSJ	10  #0099E6
  pop = EM-MUI	11 #0099E6
  pop = EMBB	12  #0099E6
  pop = EM-BOT	13  #2EBFD2
  pop = ER-h-BOT	14 #6F65FA
  pop = EXBB	15 #6F65FA
  pop = ERRB	16 #15B79D
  pop = EMQP	17  #0099E6
  pop = EMCA	18 #30099E6
  pop = EMDH	19 #3CFFF6
  pop = ER-HOO	20 #3CFFF6
  pop = ERDH	21 #3CFFF6 (02)
  pop = ERSDH	22  #3CFFF6
  pop = ERPE	23 #EBB35C
  pop = ERPA	24  #EBB35C
  pop = ERTS	25 #EBB35C
  pop = ERDM	26 #EBB35C
  pop = ERDW	27 #EBB35C
  pop = ERMK	28 #FF6600 2EBFD2  15B79D
  pop = ESTL	29 #FF6600 
  pop = ER-MOZ	30 #FF6600 
  pop = CBS	31#FF6600 
  pop = E	32 #FF6600  
  pop = PS	33 #FF6600 
  pop = SH	34 #FF6600 
  
  
  Color per pop, according to structure K=2:
    pop = EMDP	1  #0099E6
  pop = EMPN	2 #0099E6
  pop = EMKZ	3 #0099E6
  pop = EMHK	4 #0099E6
  pop = EMDB	5 #0099E6
  pop = EMJB	6 #0099E6
  pop = EMKM	7 #0099E6
  pop = EMMP	8 #0099E6
  pop = EMVM	9 #0099E6
  pop = EMSJ	10  #0099E6
  pop = EM-MUI	11 #A0099E6
  pop = EMBB	12  #40099E6
  pop = EM-BOT	13  #0099E6
  pop = ER-h-BOT	14 #0099E6
  pop = EXBB	15 #0099E6
  pop = ERRB	16 #0099E6
  pop = EMQP	17  #40099E6
  pop = EMCA	18 #0099E6
  pop = EMDH	19 #0099E6
  pop = ER-HOO	20 #0099E6
  pop = ERDH	21 #0099E6
  pop = ERSDH	22  #0099E6
  pop = ERPE	23 #FF6600
  pop = ERPA	24  #FF6600
  pop = ERTS	25 #FF6600
  pop = ERDM	26 #FF6600
  pop = ERDW	27 #FF6600
  pop = ERMK	28 #FF6600
  pop = ESTL	29 #FF6600 
  pop = ER-MOZ	30 #FF6600 
  pop = CBS	31#FF6600 
  pop = E	32 #FF6600  
  pop = PS	33 #FF6600 
  pop = SH	34 #FF6600 
  
  