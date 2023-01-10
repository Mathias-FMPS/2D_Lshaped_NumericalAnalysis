
#Path
SuiteSparsedir=/home/student/Bureau/SuiteSparse
agmgdir=/home/student/Bureau/AGMG_3.3.5-aca/SRC

# librairies de SuiteSparse
L1 = $(SuiteSparsedir)/UMFPACK/Lib/libumfpack.a
L2 = $(SuiteSparsedir)/CHOLMOD/Lib/libcholmod.a 
L3 = $(SuiteSparsedir)/AMD/Lib/libamd.a 
L4 = $(SuiteSparsedir)/CAMD/Lib/libcamd.a  
L5 = $(SuiteSparsedir)/COLAMD/Lib/libcolamd.a 
L6 = $(SuiteSparsedir)/CCOLAMD/Lib/libccolamd.a 
L7 = $(SuiteSparsedir)/metis-4.0/libmetis.a
L8 = $(SuiteSparsedir)/SuiteSparse_config/libsuitesparseconfig.a
LIB = $(L1) $(L2) $(L3) $(L4) $(L5) $(L6) $(L7) $(L8) -lm -lblas -llapack

COPT = -O3 -Wall

C=gcc

BLASLAPACK=-L/usr/lib -llapack -lblas

opt = -Ofast

export opt

# MACRO
O=object
S=src
H=header


list_c= $(agmgdir)/dagmg.o $(agmgdir)/dagmg_mumps.o
lib= -lgfortran -lm $(BLASLAPACK)

list=  $(O)/main.o $(O)/generation_m.o $(O)/verification_geometrique.o $(O)/calcul_residu.o $(O)/ecart_type.o $(O)/plotting.o $(O)/plot_comp.o  $(O)/Mesure_Rho.o $(O)/Calcul_Flux.o $(O)/Comparaison_Solveur.o $(O)/time.o $(O)/affichage.o $(O)/prob.o $(O)/modif_b.o $(O)/Rho.o $(O)/factoLU_umfpack.o $(O)/umfpack.o


default: main

clean: 
	rm $(O)/*.o 
	rm main


main: $(list) 
	cd $(agmgdir);make dseq
	$(C) $(COPT) $^ -o $@ $(list_c) $(lib) $(LIB) 

# règles spéciales liées à la récupération de la factorisation LU et à la libération mémoire via umfpack_di_free_numeric du type *void qui la contient(pour main et Comparaison_Solveur)
$(O)/main.o: main.c 
	@mkdir -p $(O)
	cc $(COPT) -c $< -o $@ -I$(SuiteSparsedir)/UMFPACK/Include \
	-I$(SuiteSparsedir)/SuiteSparse_config -I$(SuiteSparsedir)/AMD/Include -I$(H)

$(O)/Comparaison_Solveur.o: $(S)/Comparaison_Solveur.c
	cc $(COPT) -c $< -o $@ -I$(SuiteSparsedir)/UMFPACK/Include \
	-I$(SuiteSparsedir)/SuiteSparse_config -I$(SuiteSparsedir)/AMD/Include -I$(H) 

$(O)/factoLU_umfpack.o: $(S)/factoLU_umfpack.c
	cc $(COPT) -c $< -o $@ -I$(SuiteSparsedir)/UMFPACK/Include \
	-I$(SuiteSparsedir)/SuiteSparse_config -I$(SuiteSparsedir)/AMD/Include

$(O)/umfpack.o: $(S)/umfpack.c
	cc $(COPT) -c $< -o $@ -I$(SuiteSparsedir)/UMFPACK/Include \
	-I$(SuiteSparsedir)/SuiteSparse_config -I$(SuiteSparsedir)/AMD/Include

$(O)/time.o: $(S)/time.c
	cc $(COPT) -c $< -o $@


$(O)/%.o:  $(S)/%.c 
	@mkdir -p $(O)
	cc $(COPT) -c  $< -o $@  -I$(H) 


































