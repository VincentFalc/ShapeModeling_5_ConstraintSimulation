#include <iostream>
#include "gui_utils.h"
#include <igl/viewer/Viewer.h>

#include <SimulationMesh.h>
#include <FEMElementAnalytic.h>
#include <FEMElementFD.h>
#include <FEMElementAutoDiff.h>

// Viewer mesh
Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd colors_per_vertex;
Eigen::MatrixXd V_pinned;
Eigen::MatrixXd V_pinned_color;
Eigen::MatrixXd V_selected;
Eigen::MatrixXd V_selected_color;

//#define DEBUG

// sim parameters
double desiredFrameRate = 30;
bool isSimPlaying = false;
double maxStress = 0.02;

// optimization strategy
enum OptStrategy { GRADIENT_DESCENT, NEWTON };
OptStrategy optStrategy = GRADIENT_DESCENT;

// mouse interaction
enum MouseMode { NONE, PIN_NODES, APPLY_FORCE };
MouseMode mouse_mode = NONE;
int down_mouse_x = -1, down_mouse_y = -1;
int selected_v = -1;

enum ElementType { SPRING, FEM_FD, FEM_ANALYTIC, FEM_AUTODIFF };
ElementType element_type = SPRING;

// adapted from: https://stackoverflow.com/a/7811134
typedef struct { double r,g,b;} Color;
Color jetColorFromScalar(double v, double vmin, double vmax)
{
   Color c = {1.0,1.0,1.0}; // white
   double dv;

   if (v < vmin)
	  v = vmin;
   if (v > vmax)
	  v = vmax;
   dv = vmax - vmin;

   if (v < (vmin + 0.25 * dv)) {
	  c.r = 0;
	  c.g = 4 * (v - vmin) / dv;
   } else if (v < (vmin + 0.5 * dv)) {
	  c.r = 0;
	  c.b = 1 + 4 * (vmin + 0.25 * dv - v) / dv;
   } else if (v < (vmin + 0.75 * dv)) {
	  c.r = 4 * (v - vmin - 0.5 * dv) / dv;
	  c.b = 0;
   } else {
	  c.g = 1 + 4 * (vmin + 0.75 * dv - v) / dv;
	  c.b = 0;
   }

   return(c);
}

SimulationMesh simMesh;

void makeSimMesh() {
	int nRows = 30;
	int nCols = 30;
	SimulationMesh::TriangleMesh triMesh;
	SimulationMesh::generateSquareTriMesh(triMesh, -1, 0, 0.1, 0.1, nRows, nCols);
	simMesh.loadTriangleMesh(triMesh);

	if(element_type == SPRING) {
		simMesh.createSpringsFromTriangles();
	}
	else if(element_type == FEM_FD) {
		simMesh.createElementsFromTriangles<FEMElementFD>();
	}
	else if(element_type == FEM_ANALYTIC) {
		simMesh.createElementsFromTriangles<FEMElementAnalytic>();
	}
	else if(element_type == FEM_AUTODIFF ) {
		simMesh.createElementsFromTriangles<FEMElementAutoDiff>();
	}

	simMesh.addGravityForces(Vector2d(0, -9.8));
	//set some boundary conditions...
	for (int i = 0; i < nCols;i++)
	{
		Vector2d p = simMesh.x.segment<2>(2*i);
		simMesh.setPinnedNode(i, p);
	}
	for (int i = 0; i < nRows;i++)
	{
		int ii = nCols*i;
		Vector2d p = simMesh.x.segment<2>(2*ii);
		simMesh.setPinnedNode(ii, p);
	}
}

void updateViewerVertices() {
	V.resize(simMesh.numNodes, 3);
	for (int i = 0; i < simMesh.numNodes; ++i) {
		for (int j = 0; j < 2; ++j) {
			V(i,j) = simMesh.x[2*i+j];
		}
		V(i,2) = 0;
	}
}

void updateViewerVerticesColors() {
	// color vertices according to deformation
	VectorXd defEnergyPerNode;
	simMesh.computeDefoEnergyPerNode(defEnergyPerNode);
	colors_per_vertex.resize(defEnergyPerNode.size(), 3);
	for (int i = 0; i < defEnergyPerNode.size(); ++i) {
		double v = defEnergyPerNode[i];
		Color color = jetColorFromScalar(v, 0, maxStress);
		colors_per_vertex.row(i) << color.r, color.g, color.b;
	}
}

void updateGlobalPinnedVertices() {

	int nVertsPinned = simMesh.pinnedNodeElements.size();
	V_pinned.resize(nVertsPinned, 3);
	V_pinned_color.resize(nVertsPinned, 3);

	int i = 0;
	for(Element* pinnedNode : simMesh.pinnedNodeElements) {
		int nodeIndex = pinnedNode->getNodeIndex(0);
		for (int j = 0; j < 2; ++j)
			V_pinned(i,j) = simMesh.x[2*nodeIndex+j];
		V_pinned(i, 2) = 0;
		V_pinned_color.row(i) << 0, 0, 0;
		i++;
	}

}

void updateViewerFaces() {
	int nTriangles = simMesh.triangles.size();
	F.resize(nTriangles, 3);
	for (int i = 0; i < nTriangles; ++i) {
		for (int j = 0; j < 3; ++j) {
			F(i,j) = simMesh.triangles[i][2-j];
		}
	}
}

void update_display(igl::viewer::Viewer& viewer) {
	viewer.data.clear();
	viewer.data.set_mesh(V, F);
	viewer.data.compute_normals();
	viewer.data.set_colors(colors_per_vertex);

	viewer.data.set_points(V_pinned, V_pinned_color);

	if(selected_v!=-1){
		V_selected.resize(1,3); V_selected_color.resize(1,3);
		V_selected.row(0) << simMesh.x[2*selected_v+0], simMesh.x[2*selected_v+1], 0;
		V_selected_color.row(0) << 1, 0, 0;
	}
	else {
		V_selected.resize(0,0); V_selected_color.resize(0,0);
	}
	viewer.data.add_points(V_selected, V_selected_color);
}

bool callback_pre_draw(igl::viewer::Viewer& viewer) {
	if(isSimPlaying){
		if(optStrategy == GRADIENT_DESCENT)
			simMesh.solveGradientDescent();
		else if(optStrategy == NEWTON)
			simMesh.solveNewtonsMethod();
	}
	updateViewerVertices();
	updateViewerVerticesColors();
	updateGlobalPinnedVertices();
	update_display(viewer);
	return false;
}

bool callback_mouse_down(igl::viewer::Viewer& viewer, int button, int modifier) {
	down_mouse_x = viewer.current_mouse_x;
	down_mouse_y = viewer.current_mouse_y;

	if (mouse_mode == PIN_NODES) {
		int v = pick_vertex(viewer, down_mouse_x, down_mouse_y, V, F);
		if (v !=-1) {
			simMesh.togglePinnedNode(v);
			selected_v = v;
			update_display(viewer);
			return true;
		}
	}
	else if (mouse_mode == APPLY_FORCE) {
		int v = pick_vertex(viewer, down_mouse_x, down_mouse_y, V, F);
		if (v !=-1) {
			selected_v = v;
			update_display(viewer);
			return true;
		}
	}

	return false;
}

bool callback_mouse_move(igl::viewer::Viewer& viewer, int mouse_x, int mouse_y) {

	if (mouse_mode == APPLY_FORCE && selected_v != -1) {

		Eigen::Vector3d p = V.row(selected_v).transpose();
		Eigen::Vector3f pos = computePosition(viewer, mouse_x, mouse_y, V.row(selected_v));
		double k = 10;
		Eigen::Vector3d force = k*(pos.cast<double>() - p);
		simMesh.applyForceAt(selected_v, { force(0), force(1) });

		return true;
	}

	return false;
}

bool callback_mouse_up(igl::viewer::Viewer& viewer, int button, int modifier) {

	if(mouse_mode == APPLY_FORCE && selected_v != -1) {
		simMesh.applyForceAt(selected_v, {0, 0});
	}

	selected_v = -1;

	return false;
}

bool callback_init(igl::viewer::Viewer& viewer) {

	viewer.ngui->addGroup("Simulation");

	viewer.ngui->addVariable("Optimization Strategy", optStrategy, true)->setItems({ "Gradient Descent", "Newton's Method" });

	nanogui::ComboBox *widgetElement = viewer.ngui->addVariable("Element Type", element_type, true);
	widgetElement->setItems({ "spring", "FEM finite differences", "FEM analytic", "FEM AutoDiff" });
	widgetElement->setCallback([](int i){
		if(i == 0) {
			element_type = SPRING;
			simMesh.createSpringsFromTriangles();
			simMesh.x = simMesh.X;
		}
		else if(i == 1) {
			element_type = FEM_FD;
			simMesh.createElementsFromTriangles<FEMElementFD>();
			simMesh.x = simMesh.X;
		}
		else if(i == 2){
			element_type = FEM_ANALYTIC;
			simMesh.createElementsFromTriangles<FEMElementAnalytic>();
			simMesh.x = simMesh.X;
		}
		else {
			element_type = FEM_AUTODIFF;
			simMesh.createElementsFromTriangles<FEMElementAutoDiff>();
			simMesh.x = simMesh.X;
		}
	});

	// play button
	nanogui::Button *buttonPlay = viewer.ngui->addButton("Play",[]{});
	buttonPlay->setFlags(nanogui::Button::Flags::ToggleButton);
	buttonPlay->setChangeCallback([](bool pushed){isSimPlaying = pushed;});

	viewer.ngui->addVariable("Max. Stress", maxStress);

	viewer.ngui->addVariable("Mouse Mode", mouse_mode, true)->setItems({ "None","Pin Node","Apply Force" });

	viewer.ngui->addButton("Apply Force Load", []{
		simMesh.applyForceAt(1191, {-40, -40});
	});

	viewer.ngui->addButton("Test",[]{
		if(optStrategy == GRADIENT_DESCENT)
			simMesh.testGradient();
		else
			simMesh.testHessian();
	});

	viewer.screen->setSize({1200, 1000});
	viewer.screen->performLayout();
	return true;
}

int main(int argc, char *argv[]) {

	makeSimMesh();

	updateViewerVertices();
	updateViewerFaces();
	updateViewerVerticesColors();

	using Viewer = igl::viewer::Viewer;

	Viewer viewer;
	viewer.data.point_size = 10.0;
	viewer.core.align_camera_center(V,F);
	viewer.core.is_animating = true;
	viewer.callback_init = callback_init;
	viewer.callback_pre_draw = callback_pre_draw;

	viewer.callback_mouse_down = callback_mouse_down;
	viewer.callback_mouse_up = callback_mouse_up;
	viewer.callback_mouse_move = callback_mouse_move;

	viewer.launch();
}
