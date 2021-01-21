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
	// Trouve une courbe valable sur Fp avec p proche de la valeur donnée
	courbe.setRandomParam(18446744073709551615);
	// Détermine l'ordre de la courbe
	courbe.findCurveOrderHasseBSGS();
	// Trouve un point sur la courbe, qui devient le générateur du
	// sous-groupe cyclique auxquel on s'intéresse pour une utlisation cryptographique.
	courbe.findNewGen();
	// Les points d'une courbe elliptique se manipulent facilement. Ici le générateur 
	// contenu dans "courbe" ainsi que les paramètres (a, b, p) sont récupérés dans P
	// sur lequel les opérations de base s'effectuent avec aisance et simplicité
	EllPoint P(courbe.getGen());
	// On double P
	P += P;
	// On lui ajoute le générateur
	P += courbe.getGen();
	// On le multiplie
	P *= 1844674407370955;
	// On peut manipuler le générateur directement depuis la courbe. Vérifions-le :
	std::cout << (P == courbe[1844674407370955 * 3]) << std::endl; // true
	
	// Redéfinissons la courbe sur un corps plus petit et tentons de casser le DLP :
	courbe.setRandomParam(10000000, true); // 10 000 000, true : on demande à calculer l'ordre
	courbe.findNewGen();
	// Objet où le résultat est stocké (entier de taille arbitraire)
	mpz_t k; mpz_init(k);
	// On donne k l'entier où sera stocké le résultat, le point dont on veut déterminer le logarithme,
	// une limite en temps en ms, et une limite en mémoire de l'exécution de l'algorithme
	if (courbe.crackDiscreteLogBSGS(k, courbe[5000000], ~0, ~0))
		std::cout << "(DLP) Succes !" << std::endl;
	else
		std::cout << "(DLP) Echec..." << std::endl;
	
	// On peut afficher la courbe utilisée et ses paramètres ( (a, b, p), le point générateur)
	courbe.print("Exemple", true);

	// Enfin on dispose d'une fonction plot, bien pratique pour afficher une courbe par exemple
	// On fournit les tableaux des données en x, en y, les couleurs des tracés correspondants
	// et le mode de représentation : points, lignes brisées...
	// Tentons de trouver une courbe générée par au moins deux points par exemple, puis affichons
	// les sous-groupes distincts engendrés par ceux-ci :
	EllPoint G(ECParam(0,0,5,0)), H(ECParam(0,0,5,0)); // peu importe les valuers d'initialisation
	mpz_t g_order, h_order; mpz_inits(g_order, h_order, NULL);
	for (int i = 200; i < 10000; i += 10) {
		courbe.setRandomParam(i, true);
		courbe.findNewGen();
		courbe.findGenOrder();
		if (courbe.getGen().isInf() || (mpz_cmp(courbe.getGenOrder(), courbe.getCurveOrder()) == 0))
			continue; // courbe généré par un seul point ou générateur = point à l'infini (ça ne devrait pas arriver)
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
	// Si les points n'ont pas été trouvés
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