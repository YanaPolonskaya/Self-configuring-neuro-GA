
#include "GAST_KERAS.h"

using namespace std;

double The_best_NN_fitness = 0;




vector <string> map_functional = { "", "Dropout(0.1)",
"Dropout(0.2)",
"Dropout(0.3)",
"Dropout(0.4)",
"Dropout(0.5)",
"Dropout(0.6)",
"Dropout(0.7)",
"Dropout(0.8)",
"Dropout(0.9)" };



vector <string> map_activations = { "activation='softmax'",
"activation='elu'" ,
"activation='selu'",
"activation='softplus'",
"activation='softsign'" ,
"activation='relu'" ,
"activation='tanh'",
"activation='sigmoid'" ,
"activation='hard_sigmoid'",
"activation='linear'" };

GAST::GAST() //Конструктор 
{





};

std::string GAST::Get_formula(int n_ind)
{
	string formula; //Формула данного узла
	vector <double> solution;
	int n_inputs = 77;//Количество входов //?
	int contauner = 0;



	for (int j = 0; j < N; j = j + Number_of_index)
	{
		if (x2[n_ind][j + 1] != 0)
		{
			for (int k = 0; k<Number_of_index; k++)
				solution.push_back(x2[n_ind][j+k]);
				
		}
	}

	for (int i = 0; i < solution.size(); i = i + Number_of_index)
	{
		if (i == 0)
		{
			formula += "model.add(Dense(";
			contauner = 0;
			contauner = solution[i + 1];
			formula += to_string(contauner);  //?
			formula += ", input_dim=";
			formula += to_string(n_inputs);
			formula += ", ";
			contauner = 0;
			contauner = solution[i];
			formula += map_activations[contauner];
			formula += "))\n";
		}
		else
		{
			formula += "model.add(Dense(";
			contauner = 0;
			contauner = solution[i+1];
			formula += to_string(contauner);  //?
			formula += ", ";
			contauner = 0;
			contauner = solution[i];
			formula += map_activations[contauner];
			formula += "))\n";
		}

		if (solution[i + 2] != 0)
		{
			formula += "model.add(Dropout(";
			formula += to_string(solution[i + 3]);  //?
			cout << "Solutiobn i+3" << solution[i + 3]<<endl;
			formula += "))\n";
		}



	}

	
	cout << "Formula " << formula;

	return formula;


}
void GAST::Create_Python_code(std::string newfilename, int num)
{

	string PyCode;
	string oldMarkup = "\n"; //То что заменяем после Get_formula
	string Markup = "\n    "; // То на что заменяем, для правильной разметки в py коде

	ifstream file;
	ofstream fout;

	file.open("/home/ipolonskaia/projects/Linux_project_copia/Template_Code.py");
	fout.open(newfilename);
	if (file.is_open()) // вызов метода is_open()
	{
		while (!file.eof()) {

			getline(file, PyCode);

			if (PyCode.compare("    {NN_model_structure}\r")) // Находим место для вставки сгенерированного кода
				fout << PyCode;
			else //Если дошли до места, то начинаем вставлять 
			{
				PyCode = Get_formula(num);
				PyCode.insert(0, "    "); // Добавляем в начало строки пробелы для разметки
				string::size_type n = 0;
				while ((n = PyCode.find(oldMarkup, n)) != string::npos) // Находим и заменяем все oldMarkup на Markup
				{
					PyCode.replace(n, oldMarkup.size(), Markup);
					n += Markup.size();
				}
				fout << PyCode;
				//fout << "print('" + newfilename + "')\n";
			}
		}
	}
	else {

		cout << "Файл Code.py не открыт!\n\n" << endl;
	}

	file.close();
	fout.close();

}

void GAST::Calculate_fitness(int n_population) {

	string newfolder, newfile;
	char path[66];

	newfolder = "code_history";

	mkdir("code_history", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

	newfolder += "/Population_" + to_string(n_population);
	mkdir((newfolder).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);


	if (getcwd(path, sizeof(path)) == NULL)
		perror("ошибка getcwd");

	string spath = path;

	string newpath = "sys.path.append('" + spath + "/" + newfolder + "')";


	string PyCode;
	string oldMarkup = "\n"; //То что заменяем после Get_formula
	string Markup = "\n        "; // То на что заменяем, для правильной разметки в py коде

	ifstream file;
	ofstream fout;

	file.open("/home/ipolonskaia/projects/Linux_project_copia/Template_Main.py");
	string newMain = "Main_" + to_string(n_population);
	fout.open(newMain + ".py");

	if (file.is_open()) // вызов метода is_open()
	{
		while (!file.eof()) {

			getline(file, PyCode);


			if (!PyCode.compare("{newpath}\r"))
				PyCode = (newpath).c_str();
			if (PyCode.compare("        {importer}\r")) // Находим место для вставки сгенерированного кода
				fout << PyCode;
			else //Если дошли до места, то начинаем вставлять 
			{
				for (int i = 0; i < m; i++) {
					newfile = newfolder + "/Pop_" + to_string(n_population) + "_Individ_" + to_string(i) + ".py";

					Create_Python_code(newfile, i);

					newfile = "Pop_" + to_string(n_population) + "_Individ_" + to_string(i);
					PyCode = "        from " + newfile + " import MyRNN as " + "MyRNN_" + to_string(i) + "\n";
					PyCode += "        return MyRNN_" + to_string(i) + "(" + to_string(i) + ")" + "\n";
					if (i<m - 2)
						PyCode += "    elif i==" + to_string(i + 1) + ":\n";
					else if (i == m - 2)
						PyCode += "    else:\n";
					fout << PyCode;
				}


			}
		}
	}
	else {

		cout << "Файл Template_Main.py не открыт!\n\n" << endl;
	}

	file.close();
	fout.close();



	PyObject *pName, *pModule, *pDict, *pFunc, *pValue, *presult;
	// Build the name object
	pName = PyString_FromString((newMain).c_str());

	// Load the module object
	pModule = PyImport_Import(pName);

	// pDict is a borrowed reference 
	pDict = PyModule_GetDict(pModule);

	// pFunc is also a borrowed reference 
	pFunc = PyDict_GetItemString(pDict, (char*)"Calculate_fitnesses");

	if (PyCallable_Check(pFunc))
	{
		pValue = Py_BuildValue("(i)", m);
		PyErr_Print();
		presult = PyObject_CallObject(pFunc, pValue);
		PyErr_Print();
	}
	else
	{
		PyErr_Print();
	}
	//printf("Result is %d\n", PyInt_AsLong(presult));


	//PyFloat_AsDouble(presult);

	//cout << PyList_GET_SIZE(presult) << endl;

	cout << "Fitnes:   " << endl;
	//#pragma omp parallel for
	for (int i = 0; i < m; i++) {

		//newfilename = "/code_history";
		//newfilename =	"Population_" + to_string(n_population);
		//newfilename += '/';

		fi[i] = PyFloat_AsDouble((PyList_GetItem(presult, i)));
		//tree[i].Calculate(runs, newfile);
		cout << fi[i] << endl;

	}
	//#pragma omp barrier 

	// Clean up
	Py_DECREF(pModule);
	Py_DECREF(pName);
	Py_DECREF(pValue);
	/*Py_DECREF(pDict);
	Py_DECREF(pFunc);
	Py_DECREF(presult);*/

}


GAST::~GAST()
{
	if (a != 0) { delete[] a; a = 0; }
	if (b != 0) { delete[] b; b = 0; }
	if (bit != 0) { delete[] bit; bit = 0; }

	if (x1 != 0)
	{
		for (int i = 0; i < m; i++)
			delete[] x1[i];
		delete[] x1;
		x1 = 0;
	}
	if (x2 != 0)
	{
		for (int i = 0; i < m; i++)
			delete[] x2[i];
		delete[] x2;
		x2 = 0;
	}
	if (x != 0)
	{
		for (int i = 0; i < m; i++)
			delete[] x[i];
		delete[] x;
		x = 0;
	}
	if (x0 != 0)
	{
		for (int i = 0; i < m; i++)
			delete[] x0[i];
		delete[] x0;
		x0 = 0;
	}
	if (x00 != 0)
	{
		for (int i = 0; i < m; i++)
			delete[] x00[i];
		delete[] x00;
		x00 = 0;
	}

	if (fi != 0) { delete[] fi; fi = 0; }
	if (pi != 0) { delete[] pi; pi = 0; }
	if (num != 0) { delete[] num; num = 0; }

	if (opt1 != 0) { delete[] opt1; opt1 = 0; }
	if (opt2 != 0) { delete[] opt2; opt2 = 0; }
	if (opt3 != 0) { delete[] opt3; opt3 = 0; }

	if (selection_ver != 0) { delete[] selection_ver; selection_ver = 0; }
	if (crossing_ver != 0) { delete[] crossing_ver; crossing_ver = 0; }
	if (mutation_ver != 0) { delete[] mutation_ver; mutation_ver = 0; }

	if (type_of_selection != 0) { delete[] type_of_selection; type_of_selection = 0; }
	if (type_of_crossing != 0) { delete[] type_of_crossing; type_of_crossing = 0; }
	if (type_of_mutation != 0) { delete[] type_of_mutation; type_of_mutation = 0; }

	if (av_fs != 0) { delete[] av_fs; av_fs = 0; }
	if (av_fc != 0) { delete[] av_fc; av_fc = 0; }
	if (av_fm != 0) { delete[] av_fm; av_fm = 0; }

};

void GAST::a_f(int *options, int size, int *s, int numb1)
{


	int i = 0;

	for (int j = 0; j<size; j++)
	{
		if (options[j] == numb1)
		{
			s[i] = j;
			i = i + 1;
		}
	}


}

int GAST::sort(double *a, int m)

{
	double max = a[0];
	for (int i = 0; i < m; i++)
		if (a[i]>max) {
			max = a[i];
		}

	int i = 0;
	for (i = 0; i < m; i++)
	{
		if (a[i] == max) break;
	}

	return i;

}


int GAST::max1(double *a, int size)
{
	double max = a[0];
	for (int i = 0; i < size; i++)
		if (a[i]>max) {
			max = a[i];
		}

	int i = 0;
	for (i = 0; i < size; i++)
	{
		if (a[i] == max) break;
	}

	return i;

}

double GAST::gener(double *pi, int m)
{
	double ran = rand() % 10000, sum = 0; //функция-генератор
	ran /= 10000.;
	for (int i = 0; i < m; i++)
	{
		sum += pi[i];
		if (sum>ran)
			return i;
	}
}

void GAST::count(int *mass, int size, int *mas_otbor, int size_maso)
{

	for (int i = 0; i<size_maso; i++)
	{
		for (int j = 0; j<size; j++)
		{
			if (mass[j] == i) { mas_otbor[i] = mas_otbor[i] + 1; }
		}

	}
}

void GAST::average_fitness()
{

	int *s1, *s2, *s3, *s4, *s5, *c1, *c2, *c3, *c4, *m1, *m2, *m3;

	s1 = new int[opt1[0]];
	s2 = new int[opt1[1]];
	s3 = new int[opt1[2]];
	s4 = new int[opt1[3]];
	s5 = new int[opt1[4]];
	c1 = new int[opt2[0]];
	c2 = new int[opt2[1]];
	c3 = new int[opt2[2]];
	c4 = new int[opt2[3]];
	m1 = new int[opt3[0]];
	m2 = new int[opt3[1]];
	m3 = new int[opt3[2]];
	/*
	int **s1,**c1,**m1;

	s1=new int *[n_s];
	for (int i=0;i< n_s;i++)
	s1[i]=new int [opt1[i]];

	c1=new int *[n_c];        // Массив хранения популяции в двоичной форме
	for (int i=0;i< n_c;i++)
	s1[i]=new int [opt2[i]];

	m1=new int *[n_m];        // Массив хранения популяции в двоичной форме
	for (int i=0;i< n_m;i++)
	s1[i]=new int [opt3[i]];

	*/
	a_f(type_of_selection, m, s5, n_s - 1);
	a_f(type_of_selection, m, s4, n_s - 2);
	a_f(type_of_selection, m, s3, n_s - 3);
	a_f(type_of_selection, m, s2, n_s - 4);
	a_f(type_of_selection, m, s1, n_s - 5);
	a_f(type_of_crossing, m, c4, n_c - 1);          //Определяем номера элементов которые были получены при помощи данного оператора
	a_f(type_of_crossing, m, c3, n_c - 2);
	a_f(type_of_crossing, m, c2, n_c - 3);
	a_f(type_of_crossing, m, c1, n_c - 4);
	a_f(type_of_mutation, m, m3, n_m - 1);
	a_f(type_of_mutation, m, m2, n_m - 2);
	a_f(type_of_mutation, m, m1, n_m - 3);


	for (int i = 0; i<n_s; i++)
	{
		for (int j = 0; j<opt1[i]; j++)

		{
			if (i == 0) { av_fs[i] = av_fs[i] + fi[s1[j]]; }
			if (i == 1) { av_fs[i] = av_fs[i] + fi[s2[j]]; }
			if (i == 2) { av_fs[i] = av_fs[i] + fi[s3[j]]; }
			if (i == 3) { av_fs[i] = av_fs[i] + fi[s4[j]]; }
			if (i == 4) { av_fs[i] = av_fs[i] + fi[s5[j]]; }
		}
		av_fs[i] = av_fs[i] / opt1[i];
		//cout<<"avfsi: "<<av_fs[i]<<endl;

	}

	for (int i = 0; i<n_c; i++)
	{
		for (int j = 0; j<opt2[i]; j++)

		{
			if (i == 0) { av_fc[i] = av_fc[i] + fi[c1[j]]; }
			if (i == 1) { av_fc[i] = av_fc[i] + fi[c2[j]]; }
			if (i == 2) { av_fc[i] = av_fc[i] + fi[c3[j]]; }
			if (i == 3) { av_fc[i] = av_fc[i] + fi[c4[j]]; }

		}
		av_fc[i] = av_fc[i] / opt2[i];
		//cout<<"avfci: "<<av_fc[i]<<endl;


	}

	for (int i = 0; i<n_m; i++)
	{
		for (int j = 0; j<opt3[i]; j++)

		{
			if (i == 0) { av_fm[i] = av_fm[i] + fi[m1[j]]; }
			if (i == 1) { av_fm[i] = av_fm[i] + fi[m2[j]]; }
			if (i == 2) { av_fm[i] = av_fm[i] + fi[m3[j]]; }


		}
		av_fm[i] = av_fm[i] / opt3[i];
		//cout<<"avfmi: "<<av_fm[i]<<endl;


	}

	delete[] s1;
	delete[] s2;
	delete[] s3;
	delete[] s4;
	delete[] s5;
	delete[] c1;
	delete[] c2;
	delete[] c3;
	delete[] c4;
	delete[] m1;
	delete[] m2;
	delete[] m3;
}


std::vector <int> GAST::Returned_Neurons_on_each_layer()
{

	return Number_of_neurons_on_each_layers_solution;

}
int GAST::Returned_Number_of_layer_solution()
{
	return Number_of_layer_solution;

}
void GAST::Inicial(int m, int p, double Epsilon, int N, int Number_of_index, int Max_neurons, int Max_layers, int A_func)
{
	this->Max_neurons = Max_neurons;
	this->Max_layers = Max_layers;
	this->A_func = A_func;
	this->N = N;
	Number_of_layer_solution = 0;
	this->m = m; //Количество индивидов
	this->p = p; //Количество поколений
	this->Number_of_index = Number_of_index;

	a = new double[N];
	b = new double[N];

	bit = new int[N];  //Количество бит по каждой переменной (исходя из заданной области поиска по данной переменной)

	Lbit = 0;


	x1 = new double *[m];
	for (int i = 0; i< m; i++)
		x1[i] = new double[N];        // Массив хранения популяции в десятичной форме

	x2 = new double *[m];
	for (int i = 0; i< m; i++)         //Массив хранения популяции в вещественной форме
		x2[i] = new double[N];

	fi = new double[m];		     //Массив хранения пригодностей индивидов


	pi = new double[m];		     //Массив хранения вероятностей просчитанных селекцией

	num = new int[m * 2];			//Массив хранения номеров выбранных индивидов

	
									// Рассчет количества бит для каждого показателя
	for (int i = 0; i<N; i=i+Number_of_index)
	{
		for (int j = 0; j < Number_of_index; j++)
			a[i + j] = 0;
		b[i] = A_func;
		b[i+1] = Max_neurons;
		b[i + 2] = 1; //Dropout Yes or NO
		b[i + 3] = 1;//Dropout probability

	}

	Epsilon = 0.1;

	int counter = 0;

	for (int i = 0; i<N; i++)
	{
		double j = 0;
		if (counter == 3)
		{
			j = (b[i] - a[i]) / (Epsilon*0.1);
			
		}
		else
			 j = (b[i] - a[i]);

		cout << "J  " << j << endl;

		if (counter==3)
			counter = 0;
		else counter = counter + 1;

		double n = 0;


		while (j>pow(2, n))
		{
			n = n + 1;
		}

		bit[i] = n;
		if (j == 1)
			bit[i] = 1;

		cout << bit[i] << endl;
	
		Lbit = Lbit + bit[i];
	}

	cout << endl << "Lbit" << "  " << Lbit << endl << endl;



	ft_extr = 0;//значение функции в точке экстремума

	x = new double *[m];        // Массив хранения популяции в двоичной форме
	for (int i = 0; i< m; i++)
		x[i] = new double[Lbit];

	x0 = new double *[m];        // промежуточный массив хранения популяции
	for (int i = 0; i< m; i++)
		x0[i] = new double[Lbit];

	x00 = new double *[m];        // промежуточный массив хранения популяции
	for (int i = 0; i< m; i++)
		x00[i] = new double[Lbit];

	Des = 0;

	sumfi = 0;

	Zap = new double[m];

	fitness_values = new double[m];
	for (int i = 0; i < m; i++)
		fitness_values[i] = 0;

	///Для самонастройки

	numer_sel = 0;
	t_r1 = 0, t_r2 = 0;
	pi_m = new double[3]; //массив вероятностей мутации
	pi_m[0] = 0.2 / Lbit;//Слабая мутация
	pi_m[1] = 1. / Lbit;//Средняя мутация
	pi_m[2] = 5. / Lbit;//Сильная мутация

	ver_mut = 0; //вероятность мутации (которая будет в зависимости выбора типа на форме)

	n_s = 5, n_c = 4, n_m = 3;  //Количество различных типов селекции, мутации, скрещивания

	opt1 = new int[n_s];
	opt2 = new int[n_c];
	opt3 = new int[n_m];

	gen = 0;

	selection_ver = new double[n_s];
	crossing_ver = new double[n_c];
	mutation_ver = new double[n_m];

	av_fs = new double[n_s];  //average fitness selection
	av_fc = new double[n_c];  //average fitness crossover
	av_fm = new double[n_m];  //average fitness mutation

	type_of_selection = new int[m];
	type_of_crossing = new int[m];
	type_of_mutation = new int[m];

	for (int j = 0; j<n_s; j++)
	{
		opt1[j] = 0;
		av_fs[j] = 0;
	}
	for (int j = 0; j<n_c; j++)
	{
		opt2[j] = 0;
		av_fc[j] = 0;
	}                         //обнуление массивов
	for (int j = 0; j<n_m; j++)
	{
		opt3[j] = 0;
		av_fm[j] = 0;
	}

	sr_num_pokolen = 0, nomer1 = 0;
	kol_reshen = 0;//количество решений
	n_kol_reshen = 0;
	///
	//setenv("PYTHONPATH", ".", 1);

//	Py_Initialize();

};

void GAST::transformation()

{
	for (int i = 0; i<m; i++)
	{
		Des = 0;

		//Перевод индивида в 10ую систему

		int pl = bit[0], mi = 0, ii = 0;

		for (int ind = 0; ind<N; ind++)
		{
			Des = 0;


			if (ind>0)
			{
				ii = ii + 1;
				pl = pl + bit[ii];
				mi = mi + bit[ii - 1];
			}


			for (int k = 0; k<(pl - mi); k++)
			{
				double step2 = 1;// двойка в степени

				for (int u = 0; u<k; u++)
				{
					if (k == 1)
					{
						step2 = 2;
					}
					else {
						if (k == 2) step2 = 2;
						step2 = step2 * 2;

					}

				}

				Des = Des + x[i][pl - 1 - k] * step2;

			}

			x1[i][ind] = Des;

		}
	}

	//перевод в вещественную систему счисления
	double *eps1;          // массив хранения e для перевода в вещественную 
	eps1 = new double[N];
	for (int i = 0; i<N; i++)
	{
		double st = bit[i];
		eps1[i] = (b[i] - a[i]) / (pow(2, st));
		//cout<<"E "<<eps1[i]<<endl;
	}

	for (int i = 0; i < m; i++)
	{
		cout << "x2  ";
		int counter = 0;

		for (int ind = 0; ind < N; ind++)
		{
			if (counter == 0 || counter == 1)
			{
				x2[i][ind] = x1[i][ind] * eps1[ind] + a[ind];
				x2[i][ind] = floor(x2[i][ind] + 0.5);
			}
			if (counter == 2)
				x2[i][ind] = x1[i][ind];
			if (counter == 3)
				x2[i][ind] = x1[i][ind] * eps1[ind] + a[ind];
			
			if (counter == 3)
				counter = 0;
			else counter = counter + 1;

			cout << x2[i][ind] << "  ";
		}
		cout << endl;
	}

	delete[] eps1;
	
}


double GAST::evaluation(int n_population)

{

	transformation();

	
	
	double sumf = 0;

	Calculate_fitness(n_population);

	for (int i = 0; i<m; i++)
		sumf = sumf + fi[i];//просчет суммы пригодностей


	for (int i = 0; i<m; i++)
		//Оценка пригодности каждого индивида  firness=1/(1+f(x))
	{

		fitness = 0;



		if (fabs(fi[i] - ft_extr) <= Epsilon)
		{

			Zap[i] = 1;
			cout << endl << "Решение:  ";
			kol_reshen = kol_reshen + 1;
			for (int j = 0; j<N; j++)
				cout << x2[i][j] << "  ";

			cout << endl;
		}

		else { Zap[i] = 0; }

	}


	cout << "||||||||||||||||||||||||СЛЕДУЮЩЕЕ ПОКОЛЕНИЕ||||||||||||||||||||||||||||||||||||||||" << endl;
	
	return sumf;

}




double GAST::population(int n_population)
{


	for (int i = 0; i<m; i++)
	{
		Des = 0;

		for (int j = 0; j<Lbit; j++)
		{
			x[i][j] = rand() % 2;   //генерация индивидов (каждый индивид представляется бинарной строкой, куда входят значения каждой координаты)
			
			//cout<<x[i][j]<<"  ";
		}
		//	cout<<endl;
	}

	transformation();   //переводим сгенерированную популяцию в вещественную форму, для вычисления пригодности 

	double sumf = 0;

	Calculate_fitness(n_population);
	//fitnessfunc();

	for (int i = 0; i<m; i++)
		sumf = sumf + fi[i];//просчет суммы пригодностей

	return sumf;

}




void GAST::fitnessfunc()
{

	/*vector <int>  yvec;
	yvec.reserve(5);//Резервируем место
	for (int i = 5; i < 13; i++)
		yvec.push_back(i);*/

	
	string formula = Get_formula(2);
	cout << "Formula " << formula;


	PyObject *pName, *pModule, *pDict, *pFunc, *presult,*mylist_Dim, *mylist_Pop, *args;

	//char name[5];
	//printf("Enter the name pls:\n");
	//cin >> name;
	
//	PyRun_SimpleString("from RNN import someFunction");
	//PyRun_SimpleString("print(someFunction('h'))");
	printf("All is Ok\n");
	// Build the name object
	pName = PyString_FromString((char *)"RNN");
	// Load the module object
	pModule = PyImport_Import(pName);


	// pDict is a borrowed reference 
	pDict = PyModule_GetDict(pModule);


	// pFunc is also a borrowed reference 
	pFunc = PyDict_GetItemString(pDict, (char*)"someFunction");
	////////
	//pArgTuple = PyTuple_New(1); // создаем кортеж для передачи картежей в функцию
	//////
	mylist_Pop = PyList_New(m);
	// Многомерный список
	
	
	for (int j = 0; j < m; j++) {
		mylist_Dim = PyList_New(N);//2- of list 
		for (size_t i = 0; i != N; ++i) {
			PyList_SET_ITEM(mylist_Dim, i, PyFloat_FromDouble(x2[j][i]));
		}
		PyList_SET_ITEM(mylist_Pop, j, mylist_Dim);
	}
	
	args = PyTuple_Pack(1, mylist_Pop);

	if (PyCallable_Check(pFunc))
	{
		presult = PyObject_CallObject(pFunc,args);
		PyErr_Print();
	}
	else
	{
		PyErr_Print();
	}
	printf("Result is %f\n", PyFloat_AsDouble(presult));
	
	
	
	//double res = 1-PyFloat_AsDouble(presult);
	double res = 0;
	vector <double> results;  //слои которые включены в сеть (незануленные) (их номера)
	results.reserve(m);//Резервируем место
	for (int i = 0; i < m; i++) {
		results[i] = PyFloat_AsDouble((PyList_GetItem(presult, i)));
	}

	cout << "Results:" << endl;
	for (int i = 0; i < m; i++)
		cout << results[i] << " ";
	cout << endl;
	
	Py_DECREF(presult);
	// Clean up
	Py_DECREF(pModule);
	Py_DECREF(pName);

	cout << "Fitnesses: " << endl;
	for (int i = 0; i < m; i++)
	{
		fi[i] = 1. / (1 + fabs(ft_extr - (1 - results[i])));
		fitness_values[i] = 0;
		fitness_values[i] = results[i];
	}
}

void GAST::selection(int numer_sel, int potom)
{
	if (numer_sel == 0)
	{
		double average_fi = sumfi / m;//Средняя пригодность популяции

		for (int i = 0; i<m; i++)
		{
			pi[i] = fi[i] / (average_fi*m);

		}//Просчет вероятностей и Занесение вероятностей в массив

		num[potom] = gener(pi, m); //генерация вероятности при помощи функции
		num[potom + 1] = gener(pi, m);

		while (num[potom] == num[potom + 1])
		{

			num[potom + 1] = gener(pi, m);         //Проверка условий того чтобы родители не совпадали

		}

	}

	if (numer_sel == 1 || numer_sel == 2 || numer_sel == 3)
	{
		int tur = 0;
		if (numer_sel == 1)
		{
			tur = 2;
		}
		if (numer_sel == 2)
		{
			tur = 5;
		}
		if (numer_sel == 3)
		{
			tur = 9;
		}

		int *turnir, numm = 0;
		double *prigmass;
		turnir = new int[tur];  //выделение массива под турнир
		prigmass = new double[tur];//выделение массива пригодностей



		for (int plus = 0; plus<2; plus++)
		{
			for (int j = 0; j<tur; j++)
			{
				numm = rand() % m;
				turnir[j] = numm;
				numm = 0;
				prigmass[j] = fi[turnir[j]];
			}
			numm = max1(prigmass, tur);
			num[potom + plus] = turnir[numm];

			while (num[potom] == num[potom + 1])
			{

				for (int j = 0; j<tur; j++)
				{
					numm = rand() % m;
					turnir[j] = numm;
					numm = 0;
					prigmass[j] = fi[turnir[j]];

				}

				numm = max1(prigmass, tur);
				num[potom + 1] = turnir[numm];        //Проверка условий того чтобы родители не совпадали

			}
		}
		delete[] prigmass;
		delete[] turnir;

	}

	if (numer_sel == 4)

	{
		int *fi2;
		fi2 = new int[m];	 //массив номеров индивидов (индивиды расствляются по возрастанию пригодностей, порядок следования их номеров хранится в данном массиве)
		double *fi1;
		fi1 = new double[m];  //массив копия массива пригодностей

		for (int i = 0; i<m; i++)
		{
			fi1[i] = fi[i];
		}

		for (int i = 0; i<m; i++)
		{
			int p = sort(fi1, m);        //расставление по возрастанию пригодностей
			fi1[p] = 0;
			fi2[m - 1 - i] = p;
		}

		for (int i = 0; i<m; i++)
		{

			pi[fi2[i]] = ((2.*(i + 1)) / (m*(m + 1)));

		}

		//генерация вероятности при помощи функции
		num[potom] = gener(pi, m);
		num[potom + 1] = gener(pi, m);

		while (num[potom] == num[potom + 1])
		{
			num[potom + 1] = gener(pi, m);         //Проверка условий того чтобы родители не совпадали
		}

		delete[] fi1;
		delete[] fi2;
	}


}


void GAST::crossing(int numer_s, int n_par, int nu_pa1, int nu_pa2)
{
	int ii = 0;

	int t_r1 = 0, t_r2 = 0;//точки разрыва

	if (numer_s == 0)
	{

		t_r1 = 1 + rand() % (Lbit - 1);

		for (int j = 0; j<t_r1; j++)
		{
			x0[n_par][j] = x[num[nu_pa1]][j];

		}

		for (int jj = t_r1; jj<Lbit; jj++)
		{
			x0[n_par][jj] = x[num[nu_pa2]][jj];
		}

	}
	if (numer_s == 1)
	{

		int tr[2];
		tr[0] = 1 + rand() % (Lbit - 1);
		do { tr[1] = 1 + rand() % (Lbit - 1); } while (tr[0] == tr[1]);

		if (tr[0]<tr[1]) {
			t_r1 = tr[0];
			t_r2 = tr[1];
		}
		else {
			t_r2 = tr[0];
			t_r1 = tr[1];
		}


		for (int j = 0; j<t_r1; j++)
		{
			x0[n_par][j] = x[num[nu_pa1]][j];

		}

		for (int jj = t_r1; jj<t_r2; jj++)
		{
			x0[n_par][jj] = x[num[nu_pa2]][jj];
		}

		for (int jjj = t_r2; jjj<Lbit; jjj++)

		{
			x0[n_par][jjj] = x[num[nu_pa1]][jjj];

		}


	}

	if (numer_s == 2)
	{
		double ravn_gen; //число генерирующееся , и с его помощью выбирается от какого родителя ген



		for (int j = 0; j<Lbit; j++)
		{
			ravn_gen = rand() % 100;
			ravn_gen /= 100;



			if (ravn_gen<0.5)
			{
				x0[n_par][j] = x[num[nu_pa1]][j];

			}

			else {
				x0[n_par][j] = x[num[nu_pa2]][j];
			}
		}

	}

	if (numer_s == 3)
	{

		for (int k = 0; k<Lbit; k++)
			x0[n_par][k] = x[num[nu_pa1]][k];

	}

	for (int k = 0; k<Lbit; k++)
		x[n_par][k] = x0[n_par][k];

}

void GAST::mutation1(double ver_mut, int potom)
{

	for (int j = 0; j<Lbit; j++)
	{
		double ver_gen = rand() % 1000;
		ver_gen /= 1000;

		{if (ver_gen <= ver_mut)
		{
			if (x[potom][j] == 0) { x[potom][j] = 1; }
			else { x[potom][j] = 0; }
		}
		}
	}

}

int GAST::Model_Accuracy(vector<double> &Included_layers)
{
	PyObject *pName, *pModule, *pDict, *pFunc, *pValue, *presult, *pArgTuple, *pYVec;

	//char name[5];
	//printf("Enter the name pls:\n");
	//cin >> name;

	//	PyRun_SimpleString("from RNN import someFunction");
	//PyRun_SimpleString("print(someFunction('h'))");
	printf("All is Ok\n");
	// Build the name object
	pName = PyString_FromString((char *)"RNN");
	// Load the module object
	pModule = PyImport_Import(pName);


	// pDict is a borrowed reference 
	pDict = PyModule_GetDict(pModule);


	// pFunc is also a borrowed reference 
	pFunc = PyDict_GetItemString(pDict, (char*)"Model");
	////////
	pArgTuple = PyTuple_New(1); // создаем кортеж для передачи картежей в функцию
								//////

	pYVec = PyTuple_New(Included_layers.size()); //создаем кортеж с размером yvec.size()
	for (int i = 0; i < Included_layers.size(); ++i) {
		pValue = PyFloat_FromDouble(Included_layers[i]); //Создаем PyFloatObject объект из yvec[i] 

		PyTuple_SetItem(pYVec, i, pValue); //заполняем кортеж, вставляем ссылку на объект  pYvec o в позиции i кортежа, на который указывает pValue
	}

	PyTuple_SetItem(pArgTuple, 0, pYVec);

	if (PyCallable_Check(pFunc))
	{
		presult = PyObject_CallObject(pFunc, pArgTuple);
		PyErr_Print();
	}
	else
	{
		PyErr_Print();
	}
	printf("Result is %f\n", PyFloat_AsDouble(presult));

	double res = 1 - PyFloat_AsDouble(presult);
	Py_DECREF(presult);

	// Clean up
	Py_DECREF(pModule);

	Py_DECREF(pName);

	printf("RESULT OF MODEL:  %f\n:", res);
	
	cout <<endl<< "Result Structure:" << endl;
	for (int j = 0; j < Included_layers.size(); j++)
		cout << Included_layers[j] << " : ";
	cout << endl;

	return res;

}

std::vector<double> GAST::self_configuration()

{
	int Number_of_layer = 0;//Общее количество слоев
	vector <double> Included_layers;  //слои которые включены в сеть (незануленные) (их номера)

	setenv("PYTHONPATH", ".", 1);

	Py_Initialize();

	PyRun_SimpleString("import sys\nsys.path.append('/home/ipolonskaia/projects/Linux_project_copia')\nprint(sys.path)");
	
	sumfi = population(0);    //создание популяции

	printf("Population was created and evaluated");



	for (int i = 0; i<p; i++) //100 поколен
	{

		/

		for (int i = 0; i<m; i++)   
			Zap[i] = 0;


		int clon = max1(fi, m);  //поиск индивида для клонирования (на выходе получаем индивида с наибольшей пригодностью)

		

		double *clon_;
		clon_ = new double[Lbit];

		for (int j = 0; j<Lbit; j++)
			clon_[j] = x[clon][j];      //занесение индивида для клонирования в маасив, для последующего включения его в популяцию


										//

		for (int j = 0; j<n_s; j++)
		{
			opt1[j] = 0;
		}
		for (int j = 0; j<n_c; j++)
		{
			opt2[j] = 0;
		}                         //обнуление массивов
		for (int j = 0; j<n_m; j++)
		{
			opt3[j] = 0;
		}

		//Пороговые вероятности
		double Pij_s = 3. / (10 * n_s), Pij_c = 3. / (10 * n_c), Pij_m = 3. / (10 * n_m);
		//

		double subtracted = 0, sub = 0;
		int max_s_prig = 0;

		if (i == 0)  //на первой итерации
		{
			for (int j = 0; j<n_s; j++)
				selection_ver[j] = 1. / n_s;

			for (int j = 0; j<n_c - 1; j++)
				crossing_ver[j] = 0.9 / (n_c - 1);

			crossing_ver[n_c - 1] = 0.1;

			for (int j = 0; j<n_m; j++)
				mutation_ver[j] = 1. / n_m;

		}
		else
		{
			for (int j = 0; j<n_s; j++)
			{
				if ((selection_ver[j]<(Pij_s + (1. / (n_s*p)))) && (selection_ver[j]>Pij_s))
				{
					subtracted += selection_ver[j] - Pij_s;

					selection_ver[j] = Pij_s;
				}

				if (selection_ver[j]>(Pij_s + (1. / (n_s*p))))
				{
					selection_ver[j] = selection_ver[j] - (1. / (n_s*p));

					subtracted += (1. / (n_s*p));
				}

			}


			max_s_prig = 0;


			max_s_prig = max1(av_fs, n_s);  //ищем оператор с наибольшей средней пригодностью


			selection_ver[max_s_prig] = selection_ver[max_s_prig] + subtracted;  // прибавляем к оператору с наибольшей пригодностью все что отнято от остальных

			sub = 0;
			subtracted = 0;

			for (int j = 0; j<n_c; j++)
			{
				if ((crossing_ver[j]<(Pij_c + (1. / (n_c*p)))) && (crossing_ver[j]>Pij_c))
				{
					subtracted += crossing_ver[j] - Pij_c;
					crossing_ver[j] = Pij_c;

				}

				if (crossing_ver[j]>(Pij_c + (1. / (n_c*p))))
				{
					crossing_ver[j] = crossing_ver[j] - (1. / (n_c*p));
					subtracted += (1. / (n_c*p));
				}

			}

			max_s_prig = 0;


			max_s_prig = max1(av_fc, n_c);  //ищем оператор с наибольшей средней пригодностью

			crossing_ver[max_s_prig] = crossing_ver[max_s_prig] + subtracted;

			sub = 0;
			subtracted = 0;


			for (int j = 0; j<n_m; j++)
			{
				if ((mutation_ver[j]<(Pij_m + (1. / (n_m*p)))) && (mutation_ver[j]>Pij_m))
				{
					sub = mutation_ver[j] - Pij_m;
					mutation_ver[j] = Pij_m;
					subtracted += mutation_ver[j] - Pij_m;
				}

				if (mutation_ver[j]>(Pij_m + (1. / (n_m*p))))
				{
					mutation_ver[j] = mutation_ver[j] - (1. / (n_m*p));
					subtracted = subtracted + (1. / (n_m*p));
				}

			}

			max_s_prig = 0;


			max_s_prig = max1(av_fm, n_m);  //ищем оператор с наибольшей средней пригодностью


			mutation_ver[max_s_prig] = mutation_ver[max_s_prig] + subtracted;


			for (int j = 0; j<n_s; j++)
			{
				av_fs[j] = 0;
			}
			for (int j = 0; j<n_c; j++)
			{
				av_fc[j] = 0;
			}                         //обнуление массивов
			for (int j = 0; j<n_m; j++)
			{
				av_fm[j] = 0;
			}


		}

		for (int j = 0; j<m; j++)
		{
			double sum = 0;
			gen = 0;

			gen = gener(selection_ver, n_s);

			type_of_selection[j] = gen;

			gen = gener(crossing_ver, n_c);

			type_of_crossing[j] = gen;

			gen = gener(mutation_ver, n_m);

			type_of_mutation[j] = gen;

		}

		////////////Генерация новых потомков
		int num_pa = 0;
		for (int j = 0; j<m; j++)

		{
		

			selection(type_of_selection[j], num_pa);

			crossing(type_of_crossing[j], j, num_pa, num_pa + 1);

			ver_mut = pi_m[type_of_mutation[j]];

			mutation1(ver_mut, j);

			num_pa = num_pa + 2;

		}


		//Клонирование 
		for (int j = 0; j<Lbit; j++)
			x[m - 1][j] = clon_[j];

		//

		sumfi = evaluation(i+1);


		if (i == p - 1)
		{

			vector<double> solution;
			int best_fi = max1(fi, m);
			for (int j = 0; j<N; j++)
				solution.push_back(x2[best_fi][j]);
			///////
			
			Included_layers.reserve(50);//Резервируем место

			int exist = 0;
			
			for (int j = 0; j < N; j=j+ Number_of_index)
			{
				if (solution[j+1]!=0 )
				{
					for (int k = 0; k < Number_of_index; k++)
					{
						Included_layers.push_back(solution[j + k]);
						
					}
				}
			}

			cout << "Inc_LAYERS" << " in final result: " << Included_layers.size() << endl;

			cout << "Result Structure:" << endl;
			for (int j = 0; j < Included_layers.size(); j++)
				cout << Included_layers[j] << " , ";
			cout << endl;
			Number_of_layer = Included_layers.size();
			///////
			cout << "The result accuracy: " << fi[m - 1] << endl;
			int res = Model_Accuracy(Included_layers);
			

		}

		count(type_of_selection, m, opt1, n_s);
		count(type_of_crossing, m, opt2, n_c);   //каких размерностей массивы нам понадобятся для хранения результатов работы каждого из операторов
		count(type_of_mutation, m, opt3, n_m);

		average_fitness();

		int proverka = 0;
		int kkk = 0;
		for (int kk = 0; kk<m; kk++)
		{
			if (Zap[kk] == 1) {
				proverka = 1;
				

				kkk = kk;
				


			}
			//out<<kk<<": "<<Zap[kk]<<endl;


		}



		if (proverka == 1)
		{
			cout << endl << "РЕШЕНИЕ НАЙДЕНО!! на поколении " << i << endl;
			n_kol_reshen = n_kol_reshen + 1;
			nomer1 = i;

			int best_fi = max1(fi, m);

			vector<double> solution;

			for (int j = 0; j<N; j++)
				solution.push_back(x2[best_fi][j]);

			Included_layers.reserve(50);//Резервируем место

			int exist = 0;

			for (int j = 0; j < N; j = j + Number_of_index)
			{
				if (solution[j + 1] != 0)
				{
					for (int k = 0; k<Number_of_index; k++)
						Included_layers.push_back(solution[j + k]);
				}
			}

			cout << "Inc_LAYERS" << " in final result: " << Included_layers.size() << endl;

			cout << "Result Structure:" << endl;
			for (int j = 0; j < Included_layers.size(); j++)
				cout << Included_layers[j] << " : ";
			cout << endl;
			Number_of_layer = Included_layers.size();
			cout << "The result accuracy: " << fi[m - 1] << endl;

			int res = Model_Accuracy(Included_layers);


			break;
		}

		delete[] clon_;

	}

	//Starting calculation of the found model
	
	Py_Finalize();
	return Included_layers;

}