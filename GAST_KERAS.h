#ifndef GAST_H
#define GAST_H
#include <iostream>
#include <math.h>
#include <time.h>
#include <python2.7/Python.h>
#include <vector>
#include <sys/stat.h>
#include <fstream>
class GAST
{


	int N; //����������� ������������
	int m; // ��� �� ���������
	int p; //���-�� ���������
	double *a, *b; //������ � ����� �������� �� ������ ����������
	int *bit; //������ � ������� �������� ���-�� ��� ��� ������ ����������
	int Lbit; //��������� ���-�� ��� �� ���� ����������
	double ee; //��������
	double **x, **x0, **x00, **x1, **x2, *fi, *pi;
	int *num;
	double *t_extr, ft_extr;
	double Des;
	double sumfi;
	double Epsilon;
	double *Zap;
	double fitness;
	double *fitness_values;
	int Number_of_layer_solution;
	std::vector <int> Number_of_neurons_on_each_layers_solution;

	//��� �������� � ��������� ����

	int nI, nO, nIO, nIOkontrl; //����������� nI-����������� �����, nO-����������� ������, nXY-����������� ������� �������
	double mm, ll, mmkontrl, llkontrl;//�������� � ������� �������� ������
	double  **In, **Out, **Inkontrl, **Outkontrl; //����� ��������� �������, ������, � ������ ���������� ��� ������ ANN �������������� (��� ��������� � �������� �������)

												  //

												  //��� �������������
	int numer_sel;
	int t_r1, t_r2;
	double *pi_m; //������ ������������ �������
	double ver_mut; //����������� ������� (������� ����� � ����������� ������ ���� �� �����)
					//�������� �������� ��� ������������ � ������ ����������� ����������
	int n_s, n_c, n_m;
	double *selection_ver, *crossing_ver, *mutation_ver, otn;//options-������ �������� ��������� ���������� ����������
	int *opt1, *opt2, *opt3, *type_of_selection, *type_of_crossing, *type_of_mutation;
	double *av_fs, *av_fc, *av_fm;//������� ����������� ������� ��������� ������������ � ������ ������
	int Max_neurons, Max_layers, A_func, Number_of_index;
	int gen;

	double sr_num_pokolen, nomer1;
	int kol_reshen;//���������� �������
	double n_kol_reshen;//���������� �������� � ������ ���� ������� �������
						///

						//���������� ��� �������� ������
	double Ex, Ef, ex, ef;
	int running; //���������� ��������
	int file;
public:
	GAST();
	~GAST();
	void Inicial(int m, int p, double Epsilon, int N, int Number_of_index, int Max_neurons, int Max_layers, int A_func);
	void count(int *mass, int size, int *mas_otbor, int size_maso);
	void a_f(int *options, int size, int *s, int num);
	void average_fitness();
	int max1(double *a, int size);
	double gener(double *pi, int m);
	int sort(double *a, int m);
	void transformation();
	double population(int n_population);
	void selection(int numer_sel, int potom);
	void crossing(int numer_s, int n_par, int nu_pa1, int nu_pa2);
	void mutation1(double ver_mut, int potom);
	double evaluation(int n_population);
	void fitnessfunc();
	int Model_Accuracy(std::vector<double> &Structura);
	//void self_configuration(int N_, int m_, int p_, double Epsilon_, int running_,int file_,void *UndTree);
	std::vector <double> self_configuration();
	int Returned_Number_of_layer_solution();
	std::vector <int> Returned_Neurons_on_each_layer();
	std::string Get_formula(int n_ind);
	void Calculate_fitness(int n_population);
	void Create_Python_code(std::string newfilename, int num);






};







#endif
