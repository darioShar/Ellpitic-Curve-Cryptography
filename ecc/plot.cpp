#include "plot.h"
#include <cmath>
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif
#ifndef M_PI_4
#define M_PI_4 0.785398163397448309616
#endif

void testVtk()
{
	// Create a table with some points in it
	vtkSmartPointer<vtkTable> table =
		vtkSmartPointer<vtkTable>::New();

	vtkSmartPointer<vtkFloatArray> arrX =
		vtkSmartPointer<vtkFloatArray>::New();
	arrX->SetName("X Axis");
	table->AddColumn(arrX);

	vtkSmartPointer<vtkFloatArray> arrC =
		vtkSmartPointer<vtkFloatArray>::New();
	arrC->SetName("Cosine");
	table->AddColumn(arrC);

	vtkSmartPointer<vtkFloatArray> arrS =
		vtkSmartPointer<vtkFloatArray>::New();
	arrS->SetName("Sine");
	table->AddColumn(arrS);

	// Fill in the table with some example values
	int numPoints = 69;
	float inc = 7.5 / (numPoints - 1);
	table->SetNumberOfRows(numPoints);
	for (int i = 0; i < numPoints; ++i)
	{
		table->SetValue(i, 0, i * inc);
		table->SetValue(i, 1, cos(i * inc));
		table->SetValue(i, 2, sin(i * inc));
	}
	
	// Set up the view
	vtkSmartPointer<vtkContextView> view =
		vtkSmartPointer<vtkContextView>::New();
	view->GetRenderer()->SetBackground(1.0, 1.0, 1.0);

	// Add multiple line plots, setting the colors etc
	vtkSmartPointer<vtkChartXY> chart =
		vtkSmartPointer<vtkChartXY>::New();
	view->GetScene()->AddItem(chart);
	vtkPlot *line = chart->AddPlot(vtkChart::LINE);
#if VTK_MAJOR_VERSION <= 5
	line->SetInput(table, 0, 1);
#else
	line->SetInputData(table, 0, 1);
#endif
	line->SetColor(0, 255, 0, 255);
	line->SetWidth(1.0);
	line = chart->AddPlot(vtkChart::LINE);
#if VTK_MAJOR_VERSION <= 5
	line->SetInput(table, 0, 2);
#else
	line->SetInputData(table, 0, 2);
#endif
	line->SetColor(255, 0, 0, 255);
	line->SetWidth(5.0);
	// For dotted line, the line type can be from 2 to 5 for different dash/dot
	// patterns (see enum in vtkPen containing DASH_LINE, value 2):
#ifndef WIN32
	line->GetPen()->SetLineType(vtkPen::DASH_LINE);
#endif
	// (ifdef-ed out on Windows because DASH_LINE does not work on Windows
	//  machines with built-in Intel HD graphics card...)

	//view->GetRenderWindow()->SetMultiSamples(0);

	// Start interactor
	view->GetInteractor()->Initialize();
	view->GetInteractor()->Start();

}

void plot(const std::vector<std::vector<float>>& x_values, const std::vector<std::vector<float>>& y_values,
	std::vector<Vtk::Color> c, ChartType&& chartType, std::string name) {
	assert(x_values.size() == y_values.size());

	// Set up the view
	vtkSmartPointer<vtkContextView> view =
		vtkSmartPointer<vtkContextView>::New();
	view->GetRenderer()->SetBackground(1.0, 1.0, 1.0);

	// Add multiple line plots, setting the colors etc
	vtkSmartPointer<vtkChartXY> chart =
		vtkSmartPointer<vtkChartXY>::New();
	chart->SetTitle(name);
	chart->SetSize(vtkRectf(150.0f, 150.0f, 1280.0f, 720.0f));
	view->GetScene()->AddItem(chart);

	for (int i = 0; i < x_values.size(); i++) {
		std::size_t size = x_values[i].size();
		// Create a table with some points in it
		vtkSmartPointer<vtkTable> table =
			vtkSmartPointer<vtkTable>::New();

		vtkSmartPointer<vtkFloatArray> arrX =
			vtkSmartPointer<vtkFloatArray>::New();
		arrX->SetName("X Axis");
		table->AddColumn(arrX);

		vtkSmartPointer<vtkFloatArray> arrY =
			vtkSmartPointer<vtkFloatArray>::New();
		arrY->SetName(name == std::string("") ? "Y Axis" : name.c_str());
		table->AddColumn(arrY);

		// Fill in the table with some example values
		table->SetNumberOfRows(size);
		for (int j = 0; j < size; j++)
		{
			table->SetValue(j, 0, x_values[i][j]);
			table->SetValue(j, 1, y_values[i][j]);
		}

		vtkPlot *line;
		switch (chartType) {
		case ChartType::CTPOINTS:
			line = chart->AddPlot(vtkChart::POINTS);
			break;
		case ChartType::CTLINE:
			line = chart->AddPlot(vtkChart::LINE);
			break;
		case ChartType::CTBAR:
			line = chart->AddPlot(vtkChart::BAR);
			break;
		case ChartType::CTSTACKED:
			line = chart->AddPlot(vtkChart::STACKED);
			break;
		case ChartType::CTBAG:
			line = chart->AddPlot(vtkChart::BAG);
			break;
		case ChartType::CTFUNCTIONALBAG:
			line = chart->AddPlot(vtkChart::FUNCTIONALBAG);
			break;
		case ChartType::CTAREA:
			line = chart->AddPlot(vtkChart::AREA);
			break;
		}
		line->SetInputData(table, 0, 1);
		/*int a = sinf((float)M_PI_2 * (float)rand() / (float)RAND_MAX) * 255.0f;
		int b = (1.0f - sinf((float)M_PI_2 * (float)rand() / (float)RAND_MAX)) * 255.0f;
		int c = sinf(((float)M_PI_2 * (float)rand() / (float)RAND_MAX) + M_PI_4) * 255.0f;
		line->SetColor(b, c, a);*/
		line->SetColor(c[i].r, c[i].g, c[i].b);
		//line->SetWidth(0.5);
		line->GetPen()->SetLineType(vtkPen::SOLID_LINE);
	}
	// Start interactor
	view->GetInteractor()->Initialize();
	view->GetInteractor()->Start();
}


/*
Plot::Plot() : num_total(0), num_drawn(0) {
	// Create a table with some points in it
	vtkSmartPointer<vtkTable> table =
		vtkSmartPointer<vtkTable>::New();
	// Set up the view
	vtkSmartPointer<vtkContextView> view =
		vtkSmartPointer<vtkContextView>::New();
	// Add multiple line plots, setting the colors etc
	vtkSmartPointer<vtkChartXY> chart =
		vtkSmartPointer<vtkChartXY>::New();

	view->GetRenderer()->SetBackground(1.0, 1.0, 1.0);
	view->GetScene()->AddItem(chart);
}

void Plot::addData(const std::vector<float>& x_values, const std::vector<float>& y_values, std::string name) {
	assert(x_values.size() == y_values.size());
	auto size = x_values.size();
	vtkSmartPointer<vtkFloatArray> arrX = 
		vtkSmartPointer<vtkFloatArray>::New();
	vtkSmartPointer<vtkFloatArray> arrY =
		vtkSmartPointer<vtkFloatArray>::New();
	arrays.push_back(arrX); arrays.push_back(arrY);

	//table->SetNumberOfRows(std::max(table->GetNumberOfRows(), vtkIdType(size)));
	//table->SetNumberOfRows(116);

	arrX->SetName("X Axis");
	//arrX->SetArray(x_values, size, 1);
	table->AddColumn(arrX);

	arrY->SetName(name == std::string("") ? "Y Axis" : name.c_str());
	//arrY->SetArray(y_values, size, 1);
	table->AddColumn(arrY);

	for (int i = 0; i < size; i++) {
		table->SetValue(i, num_total, x_values[i]);
		table->SetValue(i, num_total + 1, y_values[i]);
	}
	num_total += 2;
}

void Plot::clear() {
	for (int i = 0; i < num_total; i++) {
		table->RemoveColumn(i);
	}
	chart->ClearPlots();
	arrays.clear();
	num_total = 0;
	num_drawn = 0;
}

void Plot::draw()
{
	while (num_total > num_drawn) {
		vtkPlot *line = chart->AddPlot(vtkChart::LINE);
		line->SetInputData(table, num_drawn, num_drawn + 1);
		line->SetColor(0, 255, 0, 255);
		line->SetWidth(1.0);
		// For dotted line, the line type can be from 2 to 5 for different dash/dot
		// patterns (see enum in vtkPen containing DASH_LINE, value 2):
		line->GetPen()->SetLineType(vtkPen::SOLID_LINE);
		num_drawn += 2;
	}
}

void Plot::display() {
	// Start interactor
	view->GetInteractor()->Initialize();
	view->GetInteractor()->Start();
}
*/