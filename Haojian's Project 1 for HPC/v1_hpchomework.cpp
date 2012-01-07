#include <iostream>
#include <string>
#include <omp.h>
#include <math.h>
#include <ctime>
#include <vector>
#include <stdio.h>
#include <iomanip>
using namespace std;

//**********Utility fuction*************//

void sleep(int mseconds)
{
	clock_t goal = mseconds + clock();
	while (goal > clock());
}

void printMainMenu()
{
	cout << "================= Project for High Performance Computing ===================" << endl;
	cout << " 1. openMP_Helloworld" << endl;
	cout << " 2. openMP_Helloworld (ordered version)" << endl;
	cout << " 3. Computing pi (serial version)" << endl;
	cout << " 4. Computing pi (parallel version)" << endl;
	cout << " 5. Matrix multiplication (serial version)" << endl;
	cout << " 6. Matrix multiplication (parallel version)" << endl;
	cout << " q. quit" << endl;
	cout << "====== End =============" << endl;
}

//**********Utility fuction ends*************//

void openMP_Helloworld()
{
	cout << "Please input the number of threads" <<endl;
	int thread_num;
	cin >> thread_num;
	cout << "We have " << thread_num << " threads" << endl;
	omp_set_num_threads(thread_num);
#pragma omp parallel firstprivate(thread_num)
	{
		//std::cout<<"Hello,World! ThreadId="<< omp_get_thread_num()<<std::endl;
		printf("Hello, world. I am thread %d of a team of %d threads. \n",omp_get_thread_num(), omp_get_num_threads());
	}
	cout << "\nBack to Menu with 'M' \n" << endl;
}

void openMP_Helloworld_Ordered()
{
	cout << "Please input the number of threads" <<endl;
	int thread_num;
	cin >> thread_num;
	cout << "We have " << thread_num << " threads" << endl;
	omp_set_num_threads(thread_num);
	int cur_index = 0;
#pragma omp parallel firstprivate(thread_num), shared(cur_index)
	{
		while(omp_get_thread_num() > cur_index)
		{
			sleep(100);
		}
		printf("Hello, world. I am thread %d of a team of %d threads. \n",omp_get_thread_num(), omp_get_num_threads());
		cur_index ++;
	}
	cout << "\nBack to Menu with 'M' \n" << endl;
}

void computing_PI_serial()
{
	int n;
	cout << "Please input the n in the formula:" <<endl;
	cin >> n;
	clock_t t1 = clock();
	long double tmp_Sum = 0;
	for(int i= 0; i<= n; i++)
	{
		if(i%2)
		{
			tmp_Sum +=(long double) -4 /((long double)(2 * i +1));
		}
		else
		{
			tmp_Sum +=(long double) 4 /((long double)(2 * i +1));
		}
	}
	clock_t t2 = clock();
	printf("Time cost = %d\n", t2 - t1);
	cout << setprecision(16)<< "The value of PI is approximated to be: "<< tmp_Sum << endl;
	long double real_pi = acos((long double) -1);
	cout << "The real value of PI is: " << real_pi << endl;
	cout << "The absolute error is: " << fabs(real_pi - tmp_Sum) << endl;
	cout << "The relative error is: " << fabs(real_pi - tmp_Sum)/ real_pi << endl;
	cout << "\nBack to Menu with 'M' \n" << endl;
	//printf("The value of PI is approximated to be: %.16f \n", tmp_Sum);
	//long double real_pi = acos((long double) -1);
	//printf("The real value of PI is : %.16f \n", real_pi);
	//printf("The absolute error is: %.16f \n", fabs(real_pi - tmp_Sum));
	//printf("The relative error is: %.16f \n", fabs(real_pi - tmp_Sum)/real_pi);

}

void computing_PI_parallel()
{
	int n;
	cout << "Please input the n in the formula:" <<endl;
	cin >> n;

	string opt;
	cout << "Please input the number of threads" <<endl;
	int thread_num;
	cin >> thread_num;
	cout << "We have " << thread_num << " threads" << endl;
	omp_set_num_threads(thread_num);

	clock_t t1 = clock();
	long double tmp_Sum = 0;

#pragma omp parallel for reduction(+: tmp_Sum)
	for(int i= 0; i<= n; i++)
	{
		if(i%2)
		{
			tmp_Sum +=(long double) -4 /((long double)(2 * i +1));
		}
		else
		{
			tmp_Sum +=(long double) 4 /((long double)(2 * i +1));
		}
	}
	clock_t t2 = clock();
	printf("Time cost = %d\n", t2 - t1);
	cout << setprecision(16)<< "The value of PI is approximated to be: "<< tmp_Sum << endl;
	long double real_pi = acos((long double) -1);
	cout << "The real value of PI is: " << real_pi << endl;
	cout << "The absolute error is: " << fabs(real_pi - tmp_Sum) << endl;
	cout << "The relative error is: " << fabs(real_pi - tmp_Sum)/ real_pi << endl;
	cout << "\nBack to Menu with 'M' \n" << endl;

	/*printf("If thread_num = %d, Time cost = %d\n", thread_num, t2 - t1);
	printf("The value of PI is approximated to be: %.16f \n", tmp_Sum);
	long double real_pi = acos((long double) -1);
	printf("The real value of PI is : %.16f \n", real_pi);
	printf("The absolute error is: %.16f \n", fabs(real_pi - tmp_Sum));
	printf("The relative error is: %.16f \n", fabs(real_pi - tmp_Sum)/real_pi);
	cout << "\nBack to Menu with 'M' \n" << endl;*/
}

void matrix_Multiplication_serial()
{
	//Initialize the n * n matrices
	int n;
	cout << "Please input the n to initialize the n * n matrices:" <<endl;
	cin >> n;
	cout << "Calculating...." <<endl;
	clock_t t1 = clock();
	vector<double> matrix_dataset_A, matrix_dataset_B, result;
	for(int i = 1; i <= n ; i++)
	{
		for(int j = 1; j <= n ; j++)
		{
			matrix_dataset_A.push_back(sin((double)(i * j)));
			matrix_dataset_B.push_back(cos((double)(i * j)));
		}
	}
	clock_t t2 = clock();
	cout << "Initialiation Part: Time Cost = " << t2- t1 << endl;
	for(int i = 0; i <n ; i ++)
	{	
		for(int j = 0; j < n ; j++)
		{
			double sum = 0;
			for(int k = 0; k < n ; k++)
			{
				sum += matrix_dataset_A[i*n+k] *  matrix_dataset_B[k*n+j];
			}
			result.push_back(sum);
		}
	}		

	clock_t t3= clock();
	cout << "Calculation Part: Time Cost = " << t3 - t2 << endl;
	if(n <10)
	{
		for(int i = 0; i < n; i++)
		{
			for(int j =0; j <n ; j++)
			{
				cout << result[i*n + j] << "  ";
			}
			cout << endl;
		}
	}

	cout << "\nBack to Menu with 'M' \n" << endl;
}

void matrix_Multiplication_paralle()
{
	//Initialize the n * n matrices
	int n;
	cout << "Please input the n to initialize the n * n matrices:" <<endl;
	cin >> n;

	string opt;
	cout << "Please input the number of threads" <<endl;
	int thread_num;
	cin >> thread_num;
	cout << "We have " << thread_num << " threads" << endl;
	omp_set_num_threads(thread_num);

	cout << "Calculating...." <<endl;
	vector<double> matrix_dataset_A, matrix_dataset_B, result;
	clock_t t1 = clock();
	for(int i = 1; i <= n ; i++)
	{
		for(int j = 1; j <= n ; j++)
		{
			matrix_dataset_A.push_back(sin((double)(i * j)));
			matrix_dataset_B.push_back(cos((double)(i * j)));
		}
	}

	#pragma omp parallel for firstprivate(result),lastprivate(result)
	for(int i = 0; i <n ; i ++)
	{	
		for(int j = 0; j < n ; j++)
		{
			double sum = 0;
			for(int k = 0; k < n ; k++)
			{
				sum += matrix_dataset_A[i*n+k] *  matrix_dataset_B[k*n+j];
			}
			result.push_back(sum);
		}
	}		
	clock_t t2 = clock();
	printf("Time cost = %d\n", t2 - t1);
	if(n <10)
	{
		for(int i = 0; i < n; i++)
		{
			for(int j =0; j <n ; j++)
			{
				cout << result[i*n + j] << "  ";
			}
			cout << endl;
		}
	}
	cout << "\nBack to Menu with 'M' \n" << endl;
}
void test()
{
	int n;
	cout << "Please input the n in the formula:" <<endl;
	cin >> n;
	clock_t t1 = clock();
	long double tmp_Sum = 0;
	for(int i= 0; i<= n; i++)
	{
		if(i%2)
		{
			tmp_Sum +=(long double) -4 /((long double)(2 * i +1));
			cout << tmp_Sum << endl;
		}
		else
		{
			tmp_Sum +=(long double) 4 /((long double)(2 * i +1));
			cout << tmp_Sum << endl;
		}
	}
	clock_t t2 = clock();
	cout << "Time cost = " << t2-t1 << endl;
	cout << setprecision(16)<< "The value of PI is approximated to be "<< tmp_Sum << endl;
	//printf("Time cost = %d\n", t2 - t1);
	//printf("The value of PI is approximated to be: %.16f \n", tmp_Sum);
	long double real_pi = acos((long double) -1);
	cout << "The real value of PI is : " << real_pi << endl;
	cout << "The absolute error is : " << fabs (real_pi - tmp_Sum) << endl;
	//printf("The real value of PI is : %.16f \n", real_pi);
	//printf("The absolute error is: %.16f \n", fabs(real_pi - tmp_Sum));
	//printf("The relative error is: %.16f \n", fabs(real_pi - tmp_Sum)/real_pi);
	cout << "\nBack to Menu with 'M' \n" << endl;

	/*
	long n;
	long double sum =0;
	cout << "Please input the n to initialize the computation:" <<endl;
	cin >> n;

	cout << "Please input the number of threads" <<endl;
	int thread_num;
	cin >> thread_num;
	cout << "We have " << thread_num << " threads" << endl;
	omp_set_num_threads(thread_num);

	#pragma omp parallel for reduction(+: sum)
	for(int i=0; i <= n; i ++)
	{
		sum += i;
	}
	cout << sum <<endl;
	*/
}

	
int main()
{
	printMainMenu();
	while (1) {
		string opt;
		getline(cin, opt);
		if (opt == "q")
			break;
		else if (opt == "m" || opt == "M" || opt == "Menu" || opt == "menu")
			printMainMenu();
		else if(opt == "1")
			openMP_Helloworld();
		else if(opt == "2")
			openMP_Helloworld_Ordered();
		else if(opt == "3")
			computing_PI_serial();
		else if(opt == "4")
			computing_PI_parallel();
		else if(opt == "5")
			matrix_Multiplication_serial();
		else if(opt == "6")
			matrix_Multiplication_paralle();
		else if(opt == "test")
			test();
	}
	
}
