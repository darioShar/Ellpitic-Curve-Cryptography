#pragma once
//#define vtkRenderingCore_AUTOINIT 4(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingFreeTypeOpenGL,vtkRenderingOpenGL2)
//#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL2) 
//#define vtkRenderingCore_AUTOINIT 2(vtkRenderingOpenGL2, vtkInteractionStyle)
#define vtkRenderingCore_AUTOINIT 3(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingOpenGL2)
#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL2)
#define vtkRenderingContext2D_AUTOINIT 1(vtkRenderingContextOpenGL2)
#include <vtkAutoInit.h>
#include <vtkVersion.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>
#include <vtkSmartPointer.h>
#include <vtkChartXY.h>
#include <vtkTable.h>
#include <vtkPlot.h>
#include <vtkFloatArray.h>
#include <vtkContextView.h>
#include <vtkContextScene.h>
#include <vtkPen.h>

static enum ChartType{
	CTLINE,
	CTPOINTS,
	CTBAR,
	CTSTACKED,
	CTBAG,
	CTFUNCTIONALBAG,
	CTAREA
};

namespace Vtk {
	struct Color {
		Color(double _r, double _g, double _b) : r(_r), g(_g), b(_b) {};
		double r, g, b;
	};
}

void extern testVtk();
void extern plot(const std::vector<std::vector<float>>& x_values, const std::vector<std::vector<float>>& y_values, 
	std::vector<Vtk::Color> c, ChartType&& chartType = ChartType::CTLINE, std::string name = "");


/*
class Plot {
public :
	Plot();

	void addData(const std::vector<float>& x_values, const std::vector<float>& y_values, std::string name = "");
	void clear();
	void draw();
	void display();

private :
	vtkSmartPointer<vtkTable>					table;
	std::vector<vtkSmartPointer<vtkFloatArray>> arrays;
	vtkSmartPointer<vtkContextView>				view; 
	vtkSmartPointer<vtkChartXY>					chart;
	std::size_t									num_total;
	std::size_t									num_drawn;
};

*/