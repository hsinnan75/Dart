.KEEP_STAT:

all: eva FluxEva SJ_Eva

Compiler	= g++
FLAGS		= -D NDEBUG -O3
LIB		= -lm

eva:		GeneralEvaluation.cpp
		$(Compiler) $(FLAGS) GeneralEvaluation.cpp -o eva

FluxEva:	FluxEvaluation.cpp
		$(Compiler) $(FLAGS) FluxEvaluation.cpp -o FluxEva

SJ_Eva:		SJ_Evaluation.cpp
		$(Compiler) $(FLAGS) SJ_Evaluation.cpp -o SJ_Eva
		
clean:
		rm -f eva FluxEva SJ_Eva
