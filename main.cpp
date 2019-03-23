#include <python2.7/Python.h>
#include <cstdio>
#include <vector>
#include <stdlib.h>
#include <iostream>
#include "GAST_KERAS.h"
#include <cmath>
using namespace std;
int main()
{
    printf("hello from Linux_project!\n");
	
	setlocale(0, "");
	srand(time(NULL));
	GAST OptNN_SGA;
	



	int m = 100; //кол-во индивидов
	int p = 100; //кол-во поколений
	int Max_number_of_neurons = 10, Max_number_of_layers = 10, A_func = 10; //Максимальное количество нейронов и слоев, количество функций
	int Number_of_index = 4; //кол-во показателей сети:  тип функции, количество нейронов, есть ли Dropout после слоя, вероятность Dropout-a
	int N = Number_of_index * Max_number_of_layers; //размерность пространства
	double Epsilin = 0.009; //Точность решения

	OptNN_SGA.Inicial(m, p, 0.0009, N, Number_of_index, Max_number_of_neurons, Max_number_of_layers, A_func);
	
	vector <double> Activ_f_solution = OptNN_SGA.self_configuration();
	/*

	setenv("PYTHONPATH", ".", 1);

	Py_Initialize();

	PyObject *pName, *pModule, *pDict, *pFunc, *pValue, *presult;

	//char name[5];
	//printf("Enter the name pls:\n");
	//cin >> name;
	PyRun_SimpleString("import sys\nsys.path.append('/home/ipolonskaia/projects/Linux_project')\nprint(sys.path)");
	PyRun_SimpleString("from RNN import someFunction");
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

	//////
	if (PyCallable_Check(pFunc))
	{
		pValue = Py_BuildValue("(z)", (char*)"ABCDEFGHIJKLMNOPQRSTUVWXYZ");
		PyErr_Print();
		printf("Let give this a shot!\n");
		presult = PyObject_CallObject(pFunc, pValue);
		PyErr_Print();
	}
	else
	{
		PyErr_Print();
	}
	printf("Result is %f\n", PyFloat_AsDouble(presult));
	Py_DECREF(pValue);

	// Clean up
	Py_DECREF(pModule);
	Py_DECREF(pName);
	double res = PyFloat_AsDouble(presult);
	// Finish the Python Interpreter
	Py_Finalize();
	cout << res << ' ' << res + res << '\n';

	*/
	
	


	return 0;
}