#include <Windows.h>
#include <iostream>
#include <string>
#include "utils.h"
#include "test.h"
#include "plot.h"


int main()
{
	srand(time(0));
	LARGE_INTEGER frequency;
	QueryPerformanceFrequency(&frequency);
	LARGE_INTEGER counterBegin;
	LARGE_INTEGER counterEnd;
	// Beginning counter :
	QueryPerformanceCounter(&counterBegin);
	//setvbuf(stdout, NULL, _IOFBF, 1 << 14); 	// setting larger stdout buffer (fast output)
	/************************************************************************************************************************/
	
	//test_1();
	//test_2();
	//test_3(103);
	//test_4(12);
	//test_hasse_bsgs_order();
	//test_hasse_naive_order();

	
	EllCurve courbe;
	// Trouve une courbe valable sur Fp avec p proche de la valeur donn�e
	courbe.setRandomParam(18446744073709551615);
	// D�termine l'ordre de la courbe
	courbe.findCurveOrderHasseBSGS();
	// Trouve un point sur la courbe, qui devient le g�n�rateur du
	// sous-groupe cyclique auxquel on s'int�resse pour une utlisation cryptographique.
	courbe.findNewGen();
	// Les points d'une courbe elliptique se manipulent facilement. Ici le g�n�rateur 
	// contenu dans "courbe" ainsi que les param�tres (a, b, p) sont r�cup�r�s dans P
	// sur lequel les op�rations de base s'effectuent avec aisance et simplicit�
	EllPoint P(courbe.getGen());
	// On double P
	P += P;
	// On lui ajoute le g�n�rateur
	P += courbe.getGen();
	// On le multiplie
	P *= 1844674407370955;
	// On peut manipuler le g�n�rateur directement depuis la courbe. V�rifions-le :
	std::cout << (P == courbe[1844674407370955 * 3]) << std::endl; // true
	
	// Red�finissons la courbe sur un corps plus petit et tentons de casser le DLP :
	courbe.setRandomParam(10000000, true); // 10 000 000, true : on demande � calculer l'ordre
	courbe.findNewGen();
	// Objet o� le r�sultat est stock� (entier de taille arbitraire)
	mpz_t k; mpz_init(k);
	// On donne k l'entier o� sera stock� le r�sultat, le point dont on veut d�terminer le logarithme,
	// une limite en temps en ms, et une limite en m�moire de l'ex�cution de l'algorithme
	if (courbe.crackDiscreteLogBSGS(k, courbe[5000000], ~0, ~0))
		std::cout << "(DLP) Succes !" << std::endl;
	else
		std::cout << "(DLP) Echec..." << std::endl;
	
	// On peut afficher la courbe utilis�e et ses param�tres ( (a, b, p), le point g�n�rateur)
	courbe.print("Exemple", true);

	// Enfin on dispose d'une fonction plot, bien pratique pour afficher une courbe par exemple
	// On fournit les tableaux des donn�es en x, en y, les couleurs des trac�s correspondants
	// et le mode de repr�sentation : points, lignes bris�es...
	// Tentons de trouver une courbe g�n�r�e par au moins deux points par exemple, puis affichons
	// les sous-groupes distincts engendr�s par ceux-ci :
	EllPoint G(ECParam(0,0,5,0)), H(ECParam(0,0,5,0)); // peu importe les valuers d'initialisation
	mpz_t g_order, h_order; mpz_inits(g_order, h_order, NULL);
	for (int i = 200; i < 10000; i += 10) {
		courbe.setRandomParam(i, true);
		courbe.findNewGen();
		courbe.findGenOrder();
		if (courbe.getGen().isInf() || (mpz_cmp(courbe.getGenOrder(), courbe.getCurveOrder()) == 0))
			continue; // courbe g�n�r� par un seul point ou g�n�rateur = point � l'infini (�a ne devrait pas arriver)
		G = courbe.getGen();
		mpz_set(g_order, courbe.getGenOrder());
		courbe.findNewGen();
		courbe.findGenOrder();
		if (mpz_cmp(g_order, courbe.getGenOrder()) != 0) {
			H = courbe.getGen();
			mpz_set(h_order, courbe.getGenOrder());
			break; // c'est bon
		}
	}
	// Si les points n'ont pas �t� trouv�s
	if (mpz_sgn(h_order) == 0) {
		std::cout << "(Plot) Echec..." << std::endl;
		mpz_clears(g_order, h_order, NULL);
		return 0;
	}

	std::vector<std::vector<float>> x(2), y(2);
	x[0].resize(mpz_get_ui(g_order));
	y[0].resize(mpz_get_ui(g_order));
	x[1].resize(mpz_get_ui(h_order));
	y[1].resize(mpz_get_ui(h_order));
	std::vector<Vtk::Color> couleurs = { Vtk::Color(1.0,0.0,0.0), Vtk::Color(0.0,0.0,1.0) };

	EllPoint tmp(G);
	for (int i = 0; i < mpz_get_ui(g_order); i++) {
		x[0][i] = (mpz_get_ui(tmp.getCoord().x));
		y[0][i] = (mpz_get_ui(tmp.getCoord().y));
		tmp += G;
	}
	tmp = H;
	for (int i = 0; i < mpz_get_ui(h_order); i++) {
		x[1][i] = (mpz_get_ui(tmp.getCoord().x));
		y[1][i] = (mpz_get_ui(tmp.getCoord().y));
		tmp += H;
	}
	// On affiche nos deux points
	G.print(); H.print();
	// Enfin, on plot
	plot(x, y, couleurs, ChartType::CTPOINTS, "Exemple");

	mpz_clears(g_order, h_order, NULL);

	/************************************************************************************************************************/
	QueryPerformanceCounter(&counterEnd);
	printf("\nProgram ended successfully. Time elapsed : ");
	printf("%.18lf", (double)(counterEnd.QuadPart - counterBegin.QuadPart) / (double)frequency.QuadPart);
	printf(" s.\nEnter any key to exit\n");
	fflush(stdout);
	// end program.
	getchar();

    return 0;
}

//std::vector<float> runTimes;


/*for (auto& v : value) {
Q = curve[v];
QueryPerformanceCounter(&counterBegin);
curve.crackDiscreteLogNaive(k, Q);
QueryPerformanceCounter(&counterEnd);
double t = (double)(counterEnd.QuadPart - counterBegin.QuadPart) / (double)frequency.QuadPart;
x.push_back(v);
y.push_back(t);
//runTimes.push_back((double)(counterEnd.QuadPart - counterBegin.QuadPart) / (double)frequency.QuadPart);
}*/

//runTimes.push_back((double)(counterEnd.QuadPart - counterBegin.QuadPart) / (double)frequency.QuadPart);

//std::for_each(runTimes.begin(), runTimes.end(), [](double time) {printf("%lf\n", time); });



/*for (int j = 0; j < 10; j++) {
courbe.findNewGen();
courbe.findGenOrder();
if (mpz_cmp(g_order, courbe.getGenOrder()) != 0) {
H = courbe.getGen();
mpz_set(h_order, courbe.getGenOrder());
break;
}
}
if (mpz_sgn(h_order) > 0) // c'est bon
break; */