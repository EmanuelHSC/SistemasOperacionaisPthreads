#pragma comment(lib,"pthreadVC2.lib")
#define HAVE_STRUCT_TIMESPEC
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <tgmath.h>
#include <time.h>

/*Definição Global*/
#define THREADS 8 /*Definição total de Threads*/

#define TOTAL_ROWS 10000 /*Total de Linhas da Matriz*/
#define TOTAL_COLS 10000 /*Total de Colunas da Matriz*/

#define MACRO_ROWS 1000 /*Limite de Linhas do Macro Bloco*/
#define MACRO_COLS 1000 /*Limite de Colunas do Macro Bloco*/

#define SEED 111

#define TOTAL_BLOCKS (TOTAL_ROWS * TOTAL_COLS) / (MACRO_ROWS * MACRO_COLS) /*Total de Macros Blocos a serem criados*/

/*Variavéis e Estruturas Globais*/
typedef struct { int max_rows; int max_cols; int min_value; int max_value; } ty_param_mat; /*Estrutura p/ parâmetros da Matriz*/
typedef struct { int row_start; int row_end; int col_start; int col_end; bool process; } ty_block_mat; /*Estrutura p/ os Macroblocos*/

double time_serial, time_paral; /*Variavel p/ cálculo do tempo gasto e speedup*/

int** mat_pric; /*Matriz Principal*/
int total_primos; /*Total de primos*/

ty_block_mat* mat_block; /*Macroblocos*/
pthread_mutex_t mutex; /*Mutex para trava de seção crítica*/

void* init_thread(void* param);

/*Funções*/
void destroy_matrix(ty_param_mat param_mat) {
	if (mat_pric != NULL){
		for (int i = 0; i < param_mat.max_rows; i++) free(mat_pric[i]);
		free(mat_pric);
	}
}

int ehPrimo(int n) {
	/*Variavel que informe se n é primo*/
	int flag = 1;
	if (n > 1) {
		/*Obtendo a raiz quadrada de n*/
		int raiz_quad = sqrt(n);
		for (int i = 2; i <= raiz_quad; i++) {
			/*Caso a divisão por i seja o resto != 0 informa que n não é primo e p/ busca*/
			if (n % i == 0) { 
				flag--; 
				break;
			}
		}
	}
	else flag--;

	return flag;
}

ty_block_mat* create_macro_blocks(ty_param_mat aux_param) {
	/*Posição de cada Macrobloco*/
	int id_block = 0;

	/*Aloca a estrutura com o total de Macroblocos criados*/
	ty_block_mat* aux_blocks = malloc(TOTAL_BLOCKS * sizeof(ty_block_mat));

	if (aux_blocks != NULL) {
		for (int i = 0; i < aux_param.max_rows; i += MACRO_ROWS) {
			for (int g = 0; g < aux_param.max_cols; g += MACRO_COLS) {
				/*Linha de inicio e fim Bloco*/
				aux_blocks[id_block].row_start = i;
				aux_blocks[id_block].row_end = i + MACRO_ROWS;

				/*Coluna de incio e fim do Bloco*/
				aux_blocks[id_block].col_start = g;
				aux_blocks[id_block].col_end = g + MACRO_COLS;

				/*Informar que o Bloco ainda não foi processado*/
				aux_blocks[id_block].process = false;

				id_block++;
			}
		}
	}
	else {
		perror("Falha ao criar a Estrutura dos Macroblocos!");
		exit(-2);
	}
	return aux_blocks;
}

void execute_search_parallel(ty_param_mat aux_param) {
	mat_block = create_macro_blocks(aux_param);

	if (mat_block != NULL) {
		total_primos = 0;

		/*Guarda o inicio da busca*/
		clock_t start_paral = clock();

		/*Inicializa o Mutex p/ travar a seção critica*/
		pthread_mutex_init(&mutex, NULL);

		/*Total de Threads a serem criadas*/
		pthread_t workers[THREADS];
		for (int i = 0; i < THREADS; i++) {
			/*Cria a Thread, caso dê algum erro ao criar é feita a verificação*/
			if (pthread_create(&workers[i], NULL, init_thread, mat_block) != 0) {
				perror("Falha na criação da Thread!");
				exit(-3);
			}
		}

		for (int g = 0; g < THREADS; g++) {
			/*Escalona as Threads*/
			pthread_join(workers[g], NULL);
		}

		/*Destroi a variavel mutex*/
		pthread_mutex_destroy(&mutex);

		/*Guarda o fim da busca*/
		clock_t end_paral = clock();

		/*Cálculo do tempo total da execução*/
		time_paral = (double)(end_paral - start_paral) / CLOCKS_PER_SEC;

		printf("Total de Primos   - Paralelo: %d\n", total_primos);
		printf("Tempo de Execucao - Paralelo: %3lf segundos\n\n", time_paral);
	}
}

void execute_search_serial(ty_param_mat aux_param) {
	/*Se estiver tudo correto com a Matriz é gerado a busca*/
	if (mat_pric != NULL){
		total_primos = 0;

		/*Guarda o tempo em que iniciou a busca*/
		clock_t start_serial = clock();
		for (int a = 0; a < aux_param.max_rows; a++) {
			for (int b = 0; b < aux_param.max_cols; b++) {
				/*Caso número seja primo retorna 1 como true*/
				if (ehPrimo(mat_pric[a][b]) == 1) total_primos++;
			}
		}
		/*Guarda o tempo em que terminou a busca*/
		clock_t end_serial = clock();

		/*Cálculo do tempo total da execução*/
		time_serial = (double)(end_serial - start_serial) / CLOCKS_PER_SEC;
		printf("Total de Primos   - Serial: %d\n", total_primos);
		printf("Tempo de Execucao - Serial: %3lf segundos\n\n", time_serial);
	}
}

ty_param_mat create_param_mat() {
	/*Parâmetros p/ criação da matriz principal*/
	ty_param_mat aux_param;
	aux_param.max_rows = TOTAL_ROWS;
	aux_param.max_cols = TOTAL_COLS;
	aux_param.min_value = 0;
	aux_param.max_value = 31999;

	return aux_param;
}

int** createMatrixRandomNumbers(ty_param_mat struct_mat) {
	/*Alocação de Memória n linhas definidos*/
	int** aux_mat = malloc(struct_mat.max_rows * sizeof(int*));

	/*Caso não tenha tido erro na Alocação de Memória*/
	if (aux_mat != NULL) {
		/*Inicialização do gerador de números aleatórios*/
		srand(SEED);
		for (int a = 0; a < struct_mat.max_rows; a++) {
			/*Alocação de Memória para n colunas definidas*/
			aux_mat[a] = malloc(struct_mat.max_cols * sizeof(int));

			if(aux_mat[a] != NULL){
				for (int b = 0; b < struct_mat.max_cols; b++)
					/*Dentro da posição a b na matriz é inserido um número aletório
					entre um intervalo definido*/
					aux_mat[a][b] = rand() % (struct_mat.max_value - struct_mat.min_value + 1) + struct_mat.min_value;
			}
			else {
				perror("Falha ao criar a Matriz Principal");
				exit(-5);
			}
		}
	}
	else {
		perror("Falha ao criar a Matriz Principal!");
		exit(-1);
	}
	/*Retorna para Matriz Principal o que foi gerado*/
	return aux_mat;
}

void* init_thread(void* param) {
	/*Caso param seja NULL ñ inicia a busca com Thread*/
	if (param != NULL) {
		/*Contador auxiliar de primos*/
		int aux_cont = 0;

		/*Faz o casting do param p/ estrutura do Bloco*/
		ty_block_mat* aux_param = (ty_block_mat*)param;
		for (int i = 0; i < TOTAL_BLOCKS; i++) {
			/*Impede que outra Thread acesso o mesmo Bloco*/
			pthread_mutex_lock(&mutex);
			if (aux_param[i].process == false) {
				/*Caso o Bloco não tenha sido processado, seta a variavél que informa
				que o Bloco está processado ou em processamento*/
				aux_param[i].process = true;
				pthread_mutex_unlock(&mutex);
				for (int a = aux_param[i].row_start; a < aux_param[i].row_end; a++) {
					for (int b = aux_param[i].col_start; b < aux_param[i].col_end; b++){
						if (ehPrimo(mat_pric[a][b]) == 1) aux_cont++;
					}
				}
			}
			/*Caso o Bloco esteja processado destrava a seção critíca*/
			else pthread_mutex_unlock(&mutex);
		}
		/*Atualiza a variavel global que informa total de primos*/
		pthread_mutex_lock(&mutex);
		total_primos += aux_cont;
		pthread_mutex_unlock(&mutex);
	}
	pthread_exit(0);
}

/*Execução do Programa*/
int main() {
	/*Criação dos parametros da Matriz Principal*/
	ty_param_mat param_mat = create_param_mat();

	/*Criação da Matriz Principal*/
	mat_pric = createMatrixRandomNumbers(param_mat);

	/*Executa - Busca Serial*/
	execute_search_serial(param_mat);

	/*Executa - Busca Paralela*/
	execute_search_parallel(param_mat);

	/*Cálculo do SpeedUp*/
	printf("SpeedUp: %3lf", (time_serial / time_paral));

	/*Destroi matriz8*/
	destroy_matrix(param_mat);

	return 0;
}
