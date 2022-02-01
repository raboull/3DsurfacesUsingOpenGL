//Author: Roman Abdoullaev
//Date: October 30th, 2021
//Course: CPSC 453
//Program Description: this program attempts to satisfy the requirements for assignment 3 of the course.

#define _USE_MATH_DEFINES
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <limits>
#include <stdio.h>
#include <math.h>
#include <functional>
#include "Geometry.h"
#include "GLDebug.h"
#include "Log.h"
#include "ShaderProgram.h"
#include "Shader.h"
#include "Texture.h"
#include "Window.h"
#include "glm/glm.hpp"
#include "glm/gtc/type_ptr.hpp"
#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/vector_angle.hpp >
using namespace std;

//setup global camera parameters
glm::mat4 model = glm::mat4(1.0f);
glm::mat4 view = glm::mat4(1.0f);
glm::mat4 proj = glm::mat4(1.0f);
glm::vec3 cameraPos = glm::vec3(0.0f, 0.0f, 3.0f);
glm::vec3 cameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
glm::vec3 cameraUp = glm::vec3(0.0f, 1.0f, 0.0f);
float cameraSpeed = 0.1f; //adjust for movement speed
//setup camera rotation parameters
float cameraSensitivity = 10.0f;//adjust for camera rotation speed
float cameraXoffset = 0;
float cameraYoffset = 0;
float yaw = -90.0f;
float pitch = 0.0f;
glm::vec2 camera_down_coord = glm::vec2(0.0f, 0.0f);//initialize mouse button down coordinates to origin
glm::vec2 camera_up_coord = glm::vec2(0.0f, 0.0f);//initialize mouse button up coordinate to origin

//this struct is used to hold the state of interactive portion of this program
struct State {
	//sceneNum and subSceneNum attributes hold the scene id that should be represented on the screen
	//sceneNum = 1 && subSceneNum = 1 2D Bezier curve edit scene.
	//sceneNum = 1 && subSceneNum = 2 2D B-Spline curve edit scene.
	//sceneNum = 2 && subSceneNum = 1 3D Bezier curve scene where user can move around.
	//sceneNum = 2 && subSceneNum = 2 3D B-Spline curve scene where user can move around.
	//sceneNum = 3 3D B-Spline surface of revolution scene where user can move around.
	//sceneNum = 4 && subSceneNum = 1 Tensor Product Bezier Surface example 1 scene where user can move around.
	//sceneNum = 4 && subSceneNum = 2 Tensor Product B-Spline Surface example 2 scene where user can move around.
	int sceneNum = 1;
	int subSceneNum = 1;
	bool pointChanged = false;
	bool movingPoint = false;
	bool movingCamera = false;
	bool wireframeToggle = false;//initially draw as a solid object
	
	//overwrite the == comparator operator for the State struct
	bool operator==(State const& other) const {
		return (sceneNum == other.sceneNum &&
			subSceneNum == other.subSceneNum &&
			pointChanged == other.pointChanged &&
			movingPoint == other.movingPoint &&
			movingCamera == other.movingCamera &&
			wireframeToggle == other.wireframeToggle);
	}
	//overwrite the != comparator operator for the State struct 
	bool operator!=(State const& other) const {
		return (sceneNum != other.sceneNum ||
			subSceneNum != other.subSceneNum ||
			pointChanged != other.pointChanged ||
			movingPoint != other.movingPoint ||
			movingCamera != other.movingCamera ||
			wireframeToggle != other.wireframeToggle);
	}
};

// We gave this code in one of the tutorials, so leaving it here too
void updateGPUGeometry(GPU_Geometry &gpuGeom, CPU_Geometry const &cpuGeom) {
	gpuGeom.bind();
	gpuGeom.setVerts(cpuGeom.verts);
	gpuGeom.setCols(cpuGeom.cols);
}

//this function creates points of a Bezier curve using the De Casteljau algorithm for provided control points
vector<glm::vec3> computeBezier(vector<glm::vec3>& control_pts) {
	//Input: control_pts is a vector that contains the control points
	//Output: a vector of coordinates for a Bezier curve

	int degree = control_pts.size() - 1;//compute the degree of the Bezier curve

	vector<glm::vec3> bezier_points;//this vector will hold the computed Bezier curve coordinates

	vector<glm::vec3> p;//make a copy of control points vector and it will be overwritten
	for (int i = 0; i < control_pts.size(); i++)
		p.push_back(control_pts[i]);

	//Compute the Bezier curve coordinates using De Casteljau algorithm in the following three loops
	for (float u = 0.0; u <= 1; u += 0.01) {
		for (int i = 1; i <= degree; i++) {
			for (int j = 0; j <= degree - i; j++) {
				int next = j + 1;
				p[j] = (1.0f - u) * p[j] + u * p[next];
			}
		}
		bezier_points.push_back(p[0]);
	}
	return bezier_points;
}

//this function creates points of an open quadratic B-Spline curve for the provided control points using the subdivision algorithm
vector<glm::vec3> computeBSpline(vector<glm::vec3> control_pts) {
	//Input: control_pts is a vector that contains the user created control points
	//Output: a vector of coordinates for a B-Spline curve
	int num_of_subdiv = 4;
	int current_subdiv = 1;
	
	vector<glm::vec3> bspline_pts;//this vector will hold the computed B-Spline curve coordinates
	
	//Compute the B-Spline curve coordinates using the Chaikin subdivision algorithm
	while (current_subdiv <= num_of_subdiv) {

		int n_d = 2 + ((control_pts.size()-2)*2) + 2;//calculate the anticipated number of points in current iteration
		bspline_pts.resize(n_d);
		
		int j = 1;//resulting subdivided curve index counter
		//perform the subdivision
		for (int i = 0; i < control_pts.size() - 1; i++) {
			bspline_pts[j] = (3.f / 4.f) * control_pts[i] + (1.f / 4.f) * control_pts[i + 1];
			bspline_pts[j + 1] = (1.f / 4.f) * control_pts[i] + (3.f / 4.f) * control_pts[i + 1];
			j = j + 2;
		}
		//include the first and the last point from the control points at the beginning of current iteration
		bspline_pts[0] = control_pts[0];
		bspline_pts[bspline_pts.size() - 1] = control_pts[control_pts.size() - 1];

		//now set the control points to be the curve created in this iteration, they might be used in the next iteration
		control_pts.resize(n_d);
		for (int i = 0; i < control_pts.size(); i++)
			control_pts[i] = bspline_pts[i];

		current_subdiv = current_subdiv + 1;
	}
	return bspline_pts;
}

//this struct holds information about mouse controls
struct MouseState {
	glm::vec2 mouse_coord;
	glm::vec2 mouse_down_coord;
	glm::vec2 mouse_up_coord;
	//left mouse button states
	bool mouse_clicked_left = false;
	bool mouse_release_left = false;
	bool mouse_held_left = false;
	//right mouse button states
	bool mouse_clicked_right = false;
	bool mouse_release_right = false;
	bool mouse_held_right = false;
	bool state_changed = true;
};

//this struct holds information about keyboard controls
struct KeyboardState {
	bool delete_pressed = false;
	bool up_pressed = false;
	bool down_pressed = false;
	bool n_pressed = false;
	bool r_pressed = false;
	bool t_pressed = false;
	bool state_changed = true;
};

//callback functions that provide peripheral input for program interaction
class MyCallbacks : public CallbackInterface {
public:
	MyCallbacks(ShaderProgram& shader, int screen_width, int screen_height) :
		shader(shader),
		screen_dimensions(screen_width, screen_height)
	{}

	//this function hadles keyboard button events
	virtual void keyCallback(int key, int scancode, int action, int mods) {
		if (key == GLFW_KEY_UP && action == GLFW_PRESS) {//increase the subScene number
			if (state.subSceneNum < 2) {
				state.subSceneNum++;
			}
		}
		if (key == GLFW_KEY_DOWN && action == GLFW_PRESS) {//decrease the subScene number
			if (state.subSceneNum > 1) {
				state.subSceneNum--;
			}
		}
		if (key == GLFW_KEY_RIGHT && action == GLFW_PRESS) {//increase the scene number
			if (state.sceneNum < 4) {
				state.sceneNum++;
			}
			else if (state.sceneNum == 4) {
				state.sceneNum = 1;
			}
		}
		if (key == GLFW_KEY_LEFT && action == GLFW_PRESS) {//decrease the scene number
			if (state.sceneNum > 1) {
				state.sceneNum--;
			}
			else if (state.sceneNum == 1) {
				state.sceneNum = 4;
			}
		}
		if (key == GLFW_KEY_DELETE && action == GLFW_PRESS) {//delete selected points
			keyboardState.delete_pressed = true;
			keyboardState.state_changed = true;
		}
		if (key == GLFW_KEY_R && action == GLFW_PRESS) {//reset the window by deleting all of the currently created control points
			keyboardState.r_pressed = true;
			keyboardState.state_changed = true;
		}
		if (key == GLFW_KEY_T && action == GLFW_PRESS) {//toggle the wireframe and solid rendering mode
			keyboardState.t_pressed = true;
			keyboardState.state_changed = true;
		}
		if (key == GLFW_KEY_W && action == GLFW_PRESS) {//move forward
			cameraPos = cameraPos + (cameraSpeed * cameraFront);
			view = glm::lookAt(cameraPos, cameraPos + cameraFront, cameraUp);
		}
		if (key == GLFW_KEY_S && action == GLFW_PRESS) {//move backwards
			cameraPos = cameraPos - (cameraSpeed * cameraFront);
			view = glm::lookAt(cameraPos, cameraPos + cameraFront, cameraUp);
		}
		if (key == GLFW_KEY_A && action == GLFW_PRESS) {//move left
			cameraPos -= glm::normalize(glm::cross(cameraFront, cameraUp)) * cameraSpeed;
			view = glm::lookAt(cameraPos, cameraPos + cameraFront, cameraUp);
		}
		if (key == GLFW_KEY_D && action == GLFW_PRESS) {//move right
			cameraPos += glm::normalize(glm::cross(cameraFront, cameraUp)) * cameraSpeed;
			view = glm::lookAt(cameraPos, cameraPos + cameraFront, cameraUp);
		}
	}

	//this function handles mouse button events
	virtual void mouseButtonCallback(int button, int action, int mods) {
		if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
			mouseState.mouse_down_coord = mouseState.mouse_coord;
			mouseState.mouse_clicked_left = true;
			mouseState.mouse_held_left = true;
		}
		else if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS) {
			mouseState.mouse_down_coord = mouseState.mouse_coord;
			mouseState.mouse_clicked_right = true;
			mouseState.mouse_held_right = true;
		}
		else if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
			mouseState.mouse_up_coord = mouseState.mouse_coord;
			mouseState.mouse_release_left = true;
			mouseState.mouse_held_left = false;
		}
		else if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_RELEASE){
			mouseState.mouse_up_coord = mouseState.mouse_coord;
			mouseState.mouse_release_right = true;
			mouseState.mouse_held_left = false;
		}
		mouseState.state_changed = true;
	}

	//this function handles mouse movement events
	virtual void cursorPosCallback(double xpos, double ypos) {
		mouseState.mouse_coord = glm::vec2(xpos, ypos);
		mouseState.mouse_coord = mouseState.mouse_coord / (screen_dimensions - 1.f);//component wide operations
		mouseState.mouse_coord *= 2;
		mouseState.mouse_coord -= 1.f;
		mouseState.mouse_coord.y = -mouseState.mouse_coord.y;//invert the y-axis coordinates
		mouseState.state_changed = true;
	}

	bool mouseStateChanged() {
		return mouseState.state_changed;
	}

	void mouseStateHandled() {
		mouseState.state_changed = false;
		mouseState.mouse_clicked_left = false;
		mouseState.mouse_clicked_right = false;
	}

	//getter for our mouse state
	const MouseState& getMouseState() {
		return mouseState;
	}

	bool keyboardStateChanged() {
		return keyboardState.state_changed;
	}

	void keyboardStatehandled() {
		keyboardState.state_changed = false;
		keyboardState.up_pressed = false;
		keyboardState.down_pressed = false;
		keyboardState.n_pressed = false;
		keyboardState.delete_pressed = false;
		keyboardState.r_pressed = false;
		keyboardState.t_pressed = false;
	}

	//getter for our keyboard state
	const KeyboardState& getKeyboardState() {
		return keyboardState;
	}

	virtual void scrollCallback(double xoffset, double yoffset) {
	}

	virtual void windowSizeCallback(int width, int height) {
		// The CallbackInterface::windowSizeCallback will call glViewport for us
		CallbackInterface::windowSizeCallback(width,  height);
	}

	State getState() {//return the State object that is stored in this MyCallbacks class
		return state;
	}

	void setState(int newSubdivisionLvl, int newSceneNum) {//update the State object that is stores in this MyCallbacksClass
		state.subSceneNum = newSubdivisionLvl;
		state.sceneNum = newSceneNum;
	}

	void setMovingPoint(bool newVal) {
		state.movingPoint = newVal;
	}

	void setMovingCamera(bool newVal) {
		state.movingCamera = newVal;
	}

	void setWireframeToggle(bool newVal) {
		state.wireframeToggle = newVal;
	}

private:
	ShaderProgram& shader;
	glm::vec2 screen_dimensions;
	MouseState mouseState;
	KeyboardState keyboardState;
	State state;	
};

void moveExistingPoint(vector<glm::vec3>& control_points, std::shared_ptr<MyCallbacks>& callback_controller,
					   glm::vec3* moveThisPoint, CPU_Geometry& controlPtsCPUGeom,
					   GPU_Geometry& controlPtsGPUGeom) {

	MouseState tempMouseState = callback_controller->getMouseState();

	if (tempMouseState.mouse_held_left == true) {//update interface while left mouse button is pressed down
		*moveThisPoint = glm::vec3(tempMouseState.mouse_coord, 0.0);
		//draw the intermediate point location
		controlPtsCPUGeom.verts.clear();
		controlPtsCPUGeom.cols.clear();
		//Place control points into appropriate CPU Geometry object
		for (int i = 0; i < control_points.size(); i++) {
			controlPtsCPUGeom.verts.push_back(control_points[i]);
		}
		//add colors to the vertices
		controlPtsCPUGeom.cols.resize(controlPtsCPUGeom.verts.size(), glm::vec3{ 1.0, 0.0, 0.0 });
		updateGPUGeometry(controlPtsGPUGeom, controlPtsCPUGeom);
	}
	else if (tempMouseState.mouse_release_left == true) {//update the final location when the left mouse button is released
		callback_controller->setMovingPoint(false);
		*moveThisPoint = glm::vec3(tempMouseState.mouse_coord, 0.0);
	}
}

void selectPointOnClick(const MouseState& mouseState, vector<glm::vec3>& control_points,
						vector<glm::vec3>& selected_control_points, State& state) {
	//copy point within 0.01 distance of the right click into a selected points vector
	float distance = -1.f;
	for (auto it = control_points.begin(); it != control_points.end();) {
		glm::vec3 current_ctrl_pt = *it;
		distance = glm::distance(glm::vec3(mouseState.mouse_down_coord, 0.f), current_ctrl_pt);
		if (distance < 0.05) {
			selected_control_points.push_back(current_ctrl_pt);
		}
		++it;//increment the iterator
	}
	state.pointChanged = true;
}

void deleteSelectedPts(vector<glm::vec3>& control_points, vector<glm::vec3>& selected_control_points, State& state) {
	vector<glm::vec3> temp_control_points;//initialize a temp control point vector that will hold only the points that are not present in selected control points vector
	for (int i = 0; i < control_points.size(); i++) {
		bool matchFound = false;
		for (int j = 0; j < selected_control_points.size(); j++) {
			if (control_points[i] == selected_control_points[j]) {
				matchFound = true;
				//std::cout << "Have a match!" << std::endl;
			}
		}
		if (matchFound == false)
			temp_control_points.push_back(control_points[i]);
	}
	control_points = temp_control_points;//update the control point vectors only with points that were not present in the selected points vector
	selected_control_points.clear();//remove the control points since they are no longer present
	state.pointChanged = true;
}

void resetPts(vector<glm::vec3>& control_points, vector<glm::vec3>& selected_control_points,
			  CPU_Geometry& controlPtsCPUGeom, CPU_Geometry& selectedControlPtsCPUGeom,
			  CPU_Geometry& curveCPUGeom, CPU_Geometry& surfaceOfRevolutionPtsCPUGeom) {

	//reset all objects that are dynamic in the program
	control_points.clear();
	selected_control_points.clear();
	controlPtsCPUGeom.verts.clear();
	controlPtsCPUGeom.cols.clear();
	selectedControlPtsCPUGeom.verts.clear();
	selectedControlPtsCPUGeom.cols.clear();
	curveCPUGeom.verts.clear();
	curveCPUGeom.cols.clear();
	surfaceOfRevolutionPtsCPUGeom.verts.clear();
	surfaceOfRevolutionPtsCPUGeom.cols.clear();
}

void setupBezierEditScene(State& state, vector<glm::vec3>& control_points, CPU_Geometry& controlPtsCPUGeom,
						  GPU_Geometry& controlPtsGPUGeom, GPU_Geometry& controlLinesGPUGeom,
						  CPU_Geometry& curveCPUGeom, GPU_Geometry& curveGPUGeom,
						  vector<glm::vec3>& selected_control_points, CPU_Geometry& selectedControlPtsCPUGeom,
						  GPU_Geometry& selectedControlPtsGPUGeom) {

	//Place control points into appropriate CPU Geometry object
	for (int i = 0; i < control_points.size(); i++) {
		controlPtsCPUGeom.verts.push_back(control_points[i]);
	}

	controlPtsCPUGeom.cols.resize(controlPtsCPUGeom.verts.size(), glm::vec3{ 1.0, 0.0, 0.0 });//add colors to the vertices
	updateGPUGeometry(controlPtsGPUGeom, controlPtsCPUGeom);

	// Reset the colors to green that will be used to draw the control lines
	controlPtsCPUGeom.cols.clear();
	controlPtsCPUGeom.cols.resize(controlPtsCPUGeom.verts.size(), glm::vec3{ 0.0, 1.0, 0.0 });

	updateGPUGeometry(controlLinesGPUGeom, controlPtsCPUGeom);

	glPointSize(20.0f);

	//now compute the points of a Bezier Curve using the provided control points
	vector<glm::vec3> bezier_points = computeBezier(control_points);
	//populare the bezierGPUGeom object with calculated bezier point coordinates
	for (int i = 0; i < bezier_points.size(); i++) {
		curveCPUGeom.verts.push_back(bezier_points[i]);
	}
	//specify the color of vertices of the bezier curve
	curveCPUGeom.cols.resize(curveCPUGeom.verts.size(), glm::vec3{ 0.0, 0.0, 1.0 });

	updateGPUGeometry(curveGPUGeom, curveCPUGeom);

	//take care of the selected points
	if (selected_control_points.size() > 0) {
		//Place selected control points into appropriate CPU Geometry object
		for (int i = 0; i < selected_control_points.size(); i++) {
			selectedControlPtsCPUGeom.verts.push_back(selected_control_points[i]);
		}
		//add black color to the vertices that are selected
		selectedControlPtsCPUGeom.cols.resize(selectedControlPtsCPUGeom.verts.size(), glm::vec3{ 0.0, 0.0, 0.0 });
		updateGPUGeometry(selectedControlPtsGPUGeom, selectedControlPtsCPUGeom);
	}
}

void setupBSplineEditScene(State& state, vector<glm::vec3>& control_points, CPU_Geometry& controlPtsCPUGeom,
	GPU_Geometry& controlPtsGPUGeom, GPU_Geometry& controlLinesGPUGeom,
	CPU_Geometry& curveCPUGeom, GPU_Geometry& curveGPUGeom,
	vector<glm::vec3>& selected_control_points, CPU_Geometry& selectedControlPtsCPUGeom,
	GPU_Geometry& selectedControlPtsGPUGeom) {

	//pass the control points to generate the Bezier curve
	//std::cout << "Inside state callback for BSPLINE scene: " << state.sceneNum << ", subscene " << state.subSceneNum << std::endl;
	//Place control points into appropriate CPU Geometry object
	for (int i = 0; i < control_points.size(); i++) {
		controlPtsCPUGeom.verts.push_back(control_points[i]);
	}

	controlPtsCPUGeom.cols.resize(controlPtsCPUGeom.verts.size(), glm::vec3{ 1.0, 0.0, 0.0 });//add colors to the vertices
	updateGPUGeometry(controlPtsGPUGeom, controlPtsCPUGeom);

	// Reset the colors to green that will be used to draw the control lines
	controlPtsCPUGeom.cols.clear();
	controlPtsCPUGeom.cols.resize(controlPtsCPUGeom.verts.size(), glm::vec3{ 0.0, 1.0, 0.0 });

	updateGPUGeometry(controlLinesGPUGeom, controlPtsCPUGeom);
	glPointSize(20.0f);

	//now compute the points of a Bezier Curve using the provided control points
	vector<glm::vec3> bSpline_points = computeBSpline(control_points);
	//populare the bezierGPUGeom object with calculated bezier point coordinates
	for (int i = 0; i < bSpline_points.size(); i++) {
		curveCPUGeom.verts.push_back(bSpline_points[i]);
	}
	//specify the color of vertices of the bezier curve
	curveCPUGeom.cols.resize(curveCPUGeom.verts.size(), glm::vec3{ 0.0, 0.0, 1.0 });

	updateGPUGeometry(curveGPUGeom, curveCPUGeom);

	//take care of the selected points
	if (selected_control_points.size() > 0) {
		//Place selected control points into appropriate CPU Geometry object
		for (int i = 0; i < selected_control_points.size(); i++) {
			selectedControlPtsCPUGeom.verts.push_back(selected_control_points[i]);
		}
		//add black color to the vertices that are selected
		selectedControlPtsCPUGeom.cols.resize(selectedControlPtsCPUGeom.verts.size(), glm::vec3{ 0.0, 0.0, 0.0 });
		updateGPUGeometry(selectedControlPtsGPUGeom, selectedControlPtsCPUGeom);
	}
}

void setupSurfaceOfRevolutionScene(vector<glm::vec3>& control_points, CPU_Geometry& controlPtsCPUGeom,
									GPU_Geometry& controlPtsGPUGeom, GPU_Geometry& controlLinesGPUGeom,
									State& state, CPU_Geometry& curveCPUGeom, GPU_Geometry& curveGPUGeom,
									vector<glm::vec3>& selected_control_points, CPU_Geometry& selectedControlPtsCPUGeom,
									GPU_Geometry& selectedControlPtsGPUGeom, CPU_Geometry& surfaceOfRevolutionTriCPUGeom,
									GPU_Geometry& surfaceOfRevolutionTriGPUGeom) {
	if (control_points.size() == 1)//if we only have one control point
		setupBezierEditScene(state, control_points, controlPtsCPUGeom, controlPtsGPUGeom, controlLinesGPUGeom, curveCPUGeom, curveGPUGeom, selected_control_points, selectedControlPtsCPUGeom, selectedControlPtsGPUGeom);
	else if (control_points.size() > 1) {//we have enough points to make a surface of revolution
		//pass the control points to generate the B-Spline curve which wil be used to create the surface of revolution
		setupBSplineEditScene(state, control_points, controlPtsCPUGeom, controlPtsGPUGeom, controlLinesGPUGeom, curveCPUGeom, curveGPUGeom, selected_control_points, selectedControlPtsCPUGeom, selectedControlPtsGPUGeom);

		vector<glm::vec3> bspline_points;//declare the B-Spline points vector
		bspline_points.resize(curveCPUGeom.verts.size());
		for (int i = 0; i < curveCPUGeom.verts.size(); i++) {//copy the B-Spline points out of the CPUgeom object
			bspline_points[i] = curveCPUGeom.verts[i];
		}

		//create a vector of rotation intervals
		vector<glm::vec1> rotate_steps;
		for (float i = 0; i < 2 * M_PI + 0.1; i = i + 0.1) {//add the extra step for the conditional check to loop back around on the surface
			rotate_steps.push_back(glm::vec1(i));
		}

		//create surface of revolution with triangles
		vector<glm::vec3> surface_of_rotation_triangles;//create a vector to hold all the surface of revolution triangle verticies
		//create the surface made of triangles
		for (int i = 0; i < rotate_steps.size() - 1; i++) {
			for (int j = 0; j < bspline_points.size() - 1; j++) {
				//bottom left coordinate of the quad
				float pt_x_bl = bspline_points[j].x * cos(rotate_steps[i].x);
				float pt_y_bl = bspline_points[j].y;
				float pt_z_bl = bspline_points[j].x * sin(rotate_steps[i].x);
				glm::vec3 pt_bl = glm::vec3(pt_x_bl, pt_y_bl, pt_z_bl);

				//top left coordinate of the quad
				float pt_x_tl = bspline_points[j + 1].x * cos(rotate_steps[i].x);
				float pt_y_tl = bspline_points[j + 1].y;
				float pt_z_tl = bspline_points[j + 1].x * sin(rotate_steps[i].x);
				glm::vec3 pt_tl = glm::vec3(pt_x_tl, pt_y_tl, pt_z_tl);

				//top right coordinate of the quad
				float pt_x_tr = bspline_points[j + 1].x * cos(rotate_steps[i + 1].x);
				float pt_y_tr = bspline_points[j + 1].y;
				float pt_z_tr = bspline_points[j + 1].x * sin(rotate_steps[i + 1].x);
				glm::vec3 pt_tr = glm::vec3(pt_x_tr, pt_y_tr, pt_z_tr);

				//bottom right coordinate of the quad
				float pt_x_br = bspline_points[j].x * cos(rotate_steps[i + 1].x);
				float pt_y_br = bspline_points[j].y;
				float pt_z_br = bspline_points[j].x * sin(rotate_steps[i + 1].x);
				glm::vec3 pt_br = glm::vec3(pt_x_br, pt_y_br, pt_z_br);

				//pack the quad triangles into a vector of points
				surface_of_rotation_triangles.push_back(pt_tr);
				surface_of_rotation_triangles.push_back(pt_tl);
				surface_of_rotation_triangles.push_back(pt_bl);
				surface_of_rotation_triangles.push_back(pt_bl);
				surface_of_rotation_triangles.push_back(pt_br);
				surface_of_rotation_triangles.push_back(pt_tr);

				//pack the quad triangles into CPU geometry object
				surfaceOfRevolutionTriCPUGeom.verts.push_back(pt_tr);
				surfaceOfRevolutionTriCPUGeom.verts.push_back(pt_tl);
				surfaceOfRevolutionTriCPUGeom.verts.push_back(pt_bl);
				surfaceOfRevolutionTriCPUGeom.verts.push_back(pt_bl);
				surfaceOfRevolutionTriCPUGeom.verts.push_back(pt_br);
				surfaceOfRevolutionTriCPUGeom.verts.push_back(pt_tr);
			}
		}
		//specify the color of surface triangles
		surfaceOfRevolutionTriCPUGeom.cols.resize(surfaceOfRevolutionTriCPUGeom.verts.size(), glm::vec3{ 1.0, 0.0, 0.0 });
		updateGPUGeometry(surfaceOfRevolutionTriGPUGeom, surfaceOfRevolutionTriCPUGeom);
		glPointSize(1.0f);
	}

}

void setUpTensorProductSurfaceEx1Scene(CPU_Geometry& exSurfCntrlPtsCPUGeom, CPU_Geometry& exSurfCntrlTriCPUGeom,
									   CPU_Geometry& exSurfTriCPUGeom, CPU_Geometry& exSurfPtsCPUGeom,
									   GPU_Geometry& exSurfCntrlPtsGPUGeom, GPU_Geometry& exSurfCntrlTriGPUGeom,
									   GPU_Geometry& exSurfPtsGPUGeom, GPU_Geometry& exSurfTriGPUGeom) {
	//Hardcoded grid of control points
	vector<glm::vec3> exSurfCntrl;//initialize the control points surface
	int cntrlDim = 5;
	//create a 5 by 5 grid of various points				
	exSurfCntrl.push_back(glm::vec3(-0.6f, -0.4f, -0.6f));//0
	exSurfCntrl.push_back(glm::vec3(-0.3f, 0.4f, -0.6f));//1
	exSurfCntrl.push_back(glm::vec3(0.0f, 0.4f, -0.6f));//2
	exSurfCntrl.push_back(glm::vec3(0.3f, 0.0f, -0.6f));//3
	exSurfCntrl.push_back(glm::vec3(0.6f, 0.1f, -0.6f));//4

	exSurfCntrl.push_back(glm::vec3(-0.6f, -0.4f, -0.3f));//5
	exSurfCntrl.push_back(glm::vec3(-0.3f, 0.4f, -0.3f));//6
	exSurfCntrl.push_back(glm::vec3(0.0f, 0.4f, -0.3f));//7
	exSurfCntrl.push_back(glm::vec3(0.3f, 0.2f, -0.3f));//8
	exSurfCntrl.push_back(glm::vec3(0.6f, 0.1f, -0.3f));//9

	exSurfCntrl.push_back(glm::vec3(-0.6f, 0.0f, 0.0f));//10
	exSurfCntrl.push_back(glm::vec3(-0.3f, 0.0f, 0.0f));//11
	exSurfCntrl.push_back(glm::vec3(0.0f, 0.0f, 0.0f));//12
	exSurfCntrl.push_back(glm::vec3(0.3f, 0.0f, 0.0f));//13
	exSurfCntrl.push_back(glm::vec3(0.6f, 0.0f, 0.0f));//14

	exSurfCntrl.push_back(glm::vec3(-0.6f, 0.0f, 0.3f));//15
	exSurfCntrl.push_back(glm::vec3(-0.3f, 0.0f, 0.3f));//16
	exSurfCntrl.push_back(glm::vec3(0.0f, 0.0f, 0.3f));//17
	exSurfCntrl.push_back(glm::vec3(0.3f, 0.0f, 0.3f));//18
	exSurfCntrl.push_back(glm::vec3(0.6f, 0.0f, 0.3f));//19

	exSurfCntrl.push_back(glm::vec3(-0.6f, 0.0f, 0.6f));//20
	exSurfCntrl.push_back(glm::vec3(-0.3f, 0.0f, 0.6f));//21
	exSurfCntrl.push_back(glm::vec3(0.0f, 0.6f, 0.6f));//22
	exSurfCntrl.push_back(glm::vec3(0.3f, 0.0f, 0.6f));//23
	exSurfCntrl.push_back(glm::vec3(0.6f, -0.8f, 0.6f));//24

	//place the control points into the corresponding CPUGeom object
	for (int i = 0; i < exSurfCntrl.size(); i++) {
		exSurfCntrlPtsCPUGeom.verts.push_back(exSurfCntrl[i]);
	}

	//create triangles to represent the surface of control points
	//Do automatically
	for (int i = 0; i < cntrlDim - 1; i++) {
		for (int j = 0; j < cntrlDim - 1; j++) {
			//Create each quad with 2 triangles
			//triangle 1
			exSurfCntrlTriCPUGeom.verts.push_back(exSurfCntrl[i + cntrlDim * j]);
			exSurfCntrlTriCPUGeom.verts.push_back(exSurfCntrl[i + 1 + cntrlDim * j]);
			exSurfCntrlTriCPUGeom.verts.push_back(exSurfCntrl[i + 1 + cntrlDim + cntrlDim * j]);
			//triangle 2
			exSurfCntrlTriCPUGeom.verts.push_back(exSurfCntrl[i + cntrlDim * j]);
			exSurfCntrlTriCPUGeom.verts.push_back(exSurfCntrl[i + 1 + cntrlDim + cntrlDim * j]);
			exSurfCntrlTriCPUGeom.verts.push_back(exSurfCntrl[i + cntrlDim + cntrlDim * j]);
		}
	}
	//create bezier curves from the control points for each row and column of control points to represent the tensor product surface
	//pull out a row of surface control points
	vector<glm::vec3> currentCntrlPtsRow;
	for (int i = 0; i < cntrlDim; i++) {
		currentCntrlPtsRow.push_back(exSurfCntrl[i]);
	}
	vector<glm::vec3> currentCntrlPtsRowBezzierCurve;
	currentCntrlPtsRowBezzierCurve = computeBezier(currentCntrlPtsRow);
	int tensorDim = currentCntrlPtsRowBezzierCurve.size();

	vector<glm::vec3> exSurfTensorSurface;
	vector<glm::vec3> exSurfTensorSurfaceCols;
	vector<glm::vec3> exSurfTensorSurfaceRows;
	vector<glm::vec3> currentCntrlPtsCol;
	vector<glm::vec3> currentCntrlPtsColBezzierCurve;

	//compute a bezier curve for each column
	for (int i = 0; i < exSurfCntrl.size(); i = i + cntrlDim) {
		currentCntrlPtsCol.clear();
		for (int j = 0; j < cntrlDim; j++)
			currentCntrlPtsCol.push_back(exSurfCntrl[j + i]);
		currentCntrlPtsColBezzierCurve.clear();
		currentCntrlPtsColBezzierCurve = computeBezier(currentCntrlPtsCol);
		for (int k = 0; k < currentCntrlPtsColBezzierCurve.size(); k++)
			exSurfTensorSurfaceCols.push_back(currentCntrlPtsColBezzierCurve[k]);
	}

	//compute bezier curve for all rows using the spacing of the row bezier curves
	for (int i = 0; i < tensorDim; i++) {
		currentCntrlPtsRow.clear();
		for (int j = 0; j < cntrlDim; j++)
			currentCntrlPtsRow.push_back(exSurfTensorSurfaceCols[j * tensorDim + i]);
		currentCntrlPtsRowBezzierCurve.clear();
		currentCntrlPtsRowBezzierCurve = computeBezier(currentCntrlPtsRow);
		for (int k = 0; k < currentCntrlPtsRowBezzierCurve.size(); k++)
			exSurfTensorSurfaceRows.push_back(currentCntrlPtsRowBezzierCurve[k]);
	}
	exSurfTensorSurface = exSurfTensorSurfaceRows;//at this point the columns contain the entire surface

	//create triangles to represent the tensor product surface points
	for (int i = 0; i < tensorDim - 1; i++) {
		for (int j = 0; j < tensorDim - 1; j++) {
			//Create each quad with 2 triangles
			//triangle 1
			exSurfTriCPUGeom.verts.push_back(exSurfTensorSurface[i + tensorDim * j]);
			exSurfTriCPUGeom.verts.push_back(exSurfTensorSurface[i + 1 + tensorDim * j]);
			exSurfTriCPUGeom.verts.push_back(exSurfTensorSurface[i + 1 + tensorDim + tensorDim * j]);
			//triangle 2
			exSurfTriCPUGeom.verts.push_back(exSurfTensorSurface[i + tensorDim * j]);
			exSurfTriCPUGeom.verts.push_back(exSurfTensorSurface[i + 1 + tensorDim + tensorDim * j]);
			exSurfTriCPUGeom.verts.push_back(exSurfTensorSurface[i + tensorDim + tensorDim * j]);
		}
	}
	//place the tensor product surface points into the corresponding CPUGeom object
	for (int i = 0; i < exSurfTensorSurfaceRows.size(); i++) {
		exSurfPtsCPUGeom.verts.push_back(exSurfTensorSurfaceRows[i]);
	}

	//specify the color of vertices of the control points
	exSurfCntrlPtsCPUGeom.cols.resize(exSurfCntrlPtsCPUGeom.verts.size(), glm::vec3{ 1.0, 0.0, 0.0 });
	//specify the color of triangles of the control point surface
	exSurfCntrlTriCPUGeom.cols.resize(exSurfCntrlTriCPUGeom.verts.size(), glm::vec3{ 0.0, 1.0, 0.0 });
	//specify the color of points of the Bezier Tensor product surface
	exSurfPtsCPUGeom.cols.resize(exSurfPtsCPUGeom.verts.size(), glm::vec3{ 0.0, 0.0, 1.0 });
	//specify the color of triangles of the control point surface
	exSurfTriCPUGeom.cols.resize(exSurfTriCPUGeom.verts.size(), glm::vec3{ 0.0, 0.0, 0.0 });

	updateGPUGeometry(exSurfCntrlPtsGPUGeom, exSurfCntrlPtsCPUGeom);
	glPointSize(10.0f);
	updateGPUGeometry(exSurfCntrlTriGPUGeom, exSurfCntrlTriCPUGeom);
	updateGPUGeometry(exSurfPtsGPUGeom, exSurfPtsCPUGeom);
	updateGPUGeometry(exSurfTriGPUGeom, exSurfTriCPUGeom);
}

void setUpTensorProductSurfaceEx2Scene(CPU_Geometry& exSurfCntrlPtsCPUGeom, CPU_Geometry& exSurfCntrlTriCPUGeom,
									   CPU_Geometry& exSurfTriCPUGeom, CPU_Geometry& exSurfPtsCPUGeom,
									   GPU_Geometry& exSurfCntrlPtsGPUGeom, GPU_Geometry& exSurfCntrlTriGPUGeom,
									   GPU_Geometry& exSurfPtsGPUGeom, GPU_Geometry& exSurfTriGPUGeom) {
	//Hardcoded grid of control points
	vector<glm::vec3> exSurfCntrl;//initialize the control points surface
	int cntrlDim = 7;
	//create a 7 by 7 grid of various points
	exSurfCntrl.push_back(glm::vec3(-0.9f, 0.0f, -0.9f));//col 1
	exSurfCntrl.push_back(glm::vec3(-0.6f, 0.0f, -0.9f));
	exSurfCntrl.push_back(glm::vec3(-0.3f, 0.0f, -0.9f));
	exSurfCntrl.push_back(glm::vec3(0.0f, 0.0f, -0.9f));
	exSurfCntrl.push_back(glm::vec3(0.3f, 0.0f, -0.9f));
	exSurfCntrl.push_back(glm::vec3(0.6f, 0.0f, -0.9f));
	exSurfCntrl.push_back(glm::vec3(0.9f, 0.0f, -0.9f));

	exSurfCntrl.push_back(glm::vec3(-0.9f, 0.0f, -0.6f));//col 2
	exSurfCntrl.push_back(glm::vec3(-0.6f, -0.4f, -0.6f));
	exSurfCntrl.push_back(glm::vec3(-0.3f, 0.4f, -0.6f));
	exSurfCntrl.push_back(glm::vec3(0.0f, 0.4f, -0.6f));
	exSurfCntrl.push_back(glm::vec3(0.3f, 0.0f, -0.6f));
	exSurfCntrl.push_back(glm::vec3(0.6f, 0.1f, -0.6f));
	exSurfCntrl.push_back(glm::vec3(0.9f, 0.0f, -0.6f));

	exSurfCntrl.push_back(glm::vec3(-0.9f, 0.0f, -0.3f));//col 3
	exSurfCntrl.push_back(glm::vec3(-0.6f, -0.4f, -0.3f));
	exSurfCntrl.push_back(glm::vec3(-0.3f, 0.4f, -0.3f));
	exSurfCntrl.push_back(glm::vec3(0.0f, 0.4f, -0.3f));
	exSurfCntrl.push_back(glm::vec3(0.3f, 0.2f, -0.3f));
	exSurfCntrl.push_back(glm::vec3(0.6f, 0.1f, -0.3f));
	exSurfCntrl.push_back(glm::vec3(0.9f, 0.0f, -0.3f));

	exSurfCntrl.push_back(glm::vec3(-0.9f, 0.0f, 0.0f));//col 4
	exSurfCntrl.push_back(glm::vec3(-0.6f, 0.0f, 0.0f));
	exSurfCntrl.push_back(glm::vec3(-0.3f, 0.0f, 0.0f));
	exSurfCntrl.push_back(glm::vec3(0.0f, 0.9f, 0.0f));
	exSurfCntrl.push_back(glm::vec3(0.3f, 0.0f, 0.0f));
	exSurfCntrl.push_back(glm::vec3(0.6f, 0.0f, 0.0f));
	exSurfCntrl.push_back(glm::vec3(0.9f, 0.0f, 0.0f));

	exSurfCntrl.push_back(glm::vec3(-0.9f, 0.0f, 0.3f));//col 5
	exSurfCntrl.push_back(glm::vec3(-0.6f, 0.0f, 0.3f));
	exSurfCntrl.push_back(glm::vec3(-0.3f, 0.0f, 0.3f));
	exSurfCntrl.push_back(glm::vec3(0.0f, 0.0f, 0.3f));
	exSurfCntrl.push_back(glm::vec3(0.3f, 0.0f, 0.3f));
	exSurfCntrl.push_back(glm::vec3(0.6f, 0.0f, 0.3f));
	exSurfCntrl.push_back(glm::vec3(0.9f, 0.0f, 0.3f));

	exSurfCntrl.push_back(glm::vec3(-0.9f, 0.0f, 0.6f));//col 6
	exSurfCntrl.push_back(glm::vec3(-0.6f, 0.0f, 0.6f));
	exSurfCntrl.push_back(glm::vec3(-0.3f, 0.0f, 0.6f));
	exSurfCntrl.push_back(glm::vec3(0.0f, 0.6f, 0.6f));
	exSurfCntrl.push_back(glm::vec3(0.3f, 0.0f, 0.6f));
	exSurfCntrl.push_back(glm::vec3(0.6f, -0.8f, 0.6f));
	exSurfCntrl.push_back(glm::vec3(0.9f, 0.0f, 0.6f));

	exSurfCntrl.push_back(glm::vec3(-0.9f, 0.0f, 0.9f));//col 7
	exSurfCntrl.push_back(glm::vec3(-0.6f, 0.0f, 0.9f));
	exSurfCntrl.push_back(glm::vec3(-0.3f, 0.0f, 0.9f));
	exSurfCntrl.push_back(glm::vec3(0.0f, 0.0f, 0.9f));
	exSurfCntrl.push_back(glm::vec3(0.3f, 0.0f, 0.9f));
	exSurfCntrl.push_back(glm::vec3(0.6f, 0.0f, 0.9f));
	exSurfCntrl.push_back(glm::vec3(0.9f, 0.0f, 0.9f));

	for (int i = 0; i < exSurfCntrl.size(); i++) {
		exSurfCntrlPtsCPUGeom.verts.push_back(exSurfCntrl[i]);
	}
	//create triangles to represent the surface of control points
	for (int i = 0; i < cntrlDim - 1; i++) {
		for (int j = 0; j < cntrlDim - 1; j++) {
			//Create each quad with 2 triangles
			//triangle 1
			exSurfCntrlTriCPUGeom.verts.push_back(exSurfCntrl[i + cntrlDim * j]);
			exSurfCntrlTriCPUGeom.verts.push_back(exSurfCntrl[i + 1 + cntrlDim * j]);
			exSurfCntrlTriCPUGeom.verts.push_back(exSurfCntrl[i + 1 + cntrlDim + cntrlDim * j]);
			//triangle 2
			exSurfCntrlTriCPUGeom.verts.push_back(exSurfCntrl[i + cntrlDim * j]);
			exSurfCntrlTriCPUGeom.verts.push_back(exSurfCntrl[i + 1 + cntrlDim + cntrlDim * j]);
			exSurfCntrlTriCPUGeom.verts.push_back(exSurfCntrl[i + cntrlDim + cntrlDim * j]);
		}
	}
	//create bezier curves from the control points for each row and column of control points to represent the tensor product surface
	//pull out a row of surface control points
	vector<glm::vec3> currentCntrlPtsRow;
	for (int i = 0; i < cntrlDim; i++) {
		currentCntrlPtsRow.push_back(exSurfCntrl[i]);
	}
	vector<glm::vec3> currentCntrlPtsRowBezzierCurve;
	currentCntrlPtsRowBezzierCurve = computeBSpline(currentCntrlPtsRow);
	int tensorDim = currentCntrlPtsRowBezzierCurve.size();

	vector<glm::vec3> exSurfTensorSurface;
	vector<glm::vec3> exSurfTensorSurfaceCols;
	vector<glm::vec3> exSurfTensorSurfaceRows;
	vector<glm::vec3> currentCntrlPtsCol;
	vector<glm::vec3> currentCntrlPtsColBezzierCurve;

	//compute a bezier curve for each column
	for (int i = 0; i < exSurfCntrl.size(); i = i + cntrlDim) {
		currentCntrlPtsCol.clear();
		for (int j = 0; j < cntrlDim; j++)
			currentCntrlPtsCol.push_back(exSurfCntrl[j + i]);
		currentCntrlPtsColBezzierCurve.clear();
		currentCntrlPtsColBezzierCurve = computeBSpline(currentCntrlPtsCol);
		for (int k = 0; k < currentCntrlPtsColBezzierCurve.size(); k++)
			exSurfTensorSurfaceCols.push_back(currentCntrlPtsColBezzierCurve[k]);
	}
	//compute bezier curve for all rows using the spacing of the row bezier curves
	for (int i = 0; i < tensorDim; i++) {
		currentCntrlPtsRow.clear();
		for (int j = 0; j < cntrlDim; j++)
			currentCntrlPtsRow.push_back(exSurfTensorSurfaceCols[j * tensorDim + i]);
		currentCntrlPtsRowBezzierCurve.clear();
		currentCntrlPtsRowBezzierCurve = computeBSpline(currentCntrlPtsRow);
		for (int k = 0; k < currentCntrlPtsRowBezzierCurve.size(); k++)
			exSurfTensorSurfaceRows.push_back(currentCntrlPtsRowBezzierCurve[k]);
	}
	exSurfTensorSurface = exSurfTensorSurfaceRows;//at this point the columns contain the entire surface

	//create triangles to represent the tensor product surface points
	for (int i = 0; i < tensorDim - 1; i++) {
		for (int j = 0; j < tensorDim - 1; j++) {
			//Create each quad with 2 triangles
			//triangle 1
			exSurfTriCPUGeom.verts.push_back(exSurfTensorSurface[i + tensorDim * j]);
			exSurfTriCPUGeom.verts.push_back(exSurfTensorSurface[i + 1 + tensorDim * j]);
			exSurfTriCPUGeom.verts.push_back(exSurfTensorSurface[i + 1 + tensorDim + tensorDim * j]);
			//triangle 2
			exSurfTriCPUGeom.verts.push_back(exSurfTensorSurface[i + tensorDim * j]);
			exSurfTriCPUGeom.verts.push_back(exSurfTensorSurface[i + 1 + tensorDim + tensorDim * j]);
			exSurfTriCPUGeom.verts.push_back(exSurfTensorSurface[i + tensorDim + tensorDim * j]);
		}
	}

	//place the tensor product surface points into the corresponding CPUGeom object
	for (int i = 0; i < exSurfTensorSurfaceRows.size(); i++) {
		exSurfPtsCPUGeom.verts.push_back(exSurfTensorSurfaceRows[i]);
	}

	//specify the color of vertices of the control points
	exSurfCntrlPtsCPUGeom.cols.resize(exSurfCntrlPtsCPUGeom.verts.size(), glm::vec3{ 1.0, 0.0, 0.0 });
	//specify the color of triangles of the control point surface
	exSurfCntrlTriCPUGeom.cols.resize(exSurfCntrlTriCPUGeom.verts.size(), glm::vec3{ 0.0, 1.0, 0.0 });
	//specify the color of points of the Bezier Tensor product surface
	exSurfPtsCPUGeom.cols.resize(exSurfPtsCPUGeom.verts.size(), glm::vec3{ 0.0, 0.0, 1.0 });
	//specify the color of triangles of the control point surface
	exSurfTriCPUGeom.cols.resize(exSurfTriCPUGeom.verts.size(), glm::vec3{ 0.0, 0.0, 0.0 });

	updateGPUGeometry(exSurfCntrlPtsGPUGeom, exSurfCntrlPtsCPUGeom);
	glPointSize(10.0f);
	updateGPUGeometry(exSurfCntrlTriGPUGeom, exSurfCntrlTriCPUGeom);
	updateGPUGeometry(exSurfPtsGPUGeom, exSurfPtsCPUGeom);
	updateGPUGeometry(exSurfTriGPUGeom, exSurfTriCPUGeom);
}

void toggleWireframeSolidView(State& state, std::shared_ptr<MyCallbacks>& callback_controller) {
	if (state.wireframeToggle == true) {
		callback_controller->setWireframeToggle(false);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}
	else if (state.wireframeToggle == false) {
		callback_controller->setWireframeToggle(true);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	}
}

int main() {
	Log::debug("Starting main");

	// WINDOW
	glfwInit();
	int screen_width = 800;
	int screen_height = 800;
	Window window(screen_width, screen_height, "CPSC 453 - A3"); // can set callbacks at construction if desired

	GLDebug::enable();

	//SHADERS
	ShaderProgram shader("shaders/test.vert", "shaders/test.frag");

	//Create and initialize our state object
	State state;
	state.sceneNum = 0;//initilize to a non standard value to trigger the draw of scene 1 subscene 1 at initialization
	state.subSceneNum = 1;

	// CALLBACKS
	std::shared_ptr<MyCallbacks> callback_controller = std::make_shared<MyCallbacks>(shader, screen_width, screen_height);
	window.setCallbacks(callback_controller); // can also update callbacks to new ones
	
	// The current CPU_Geometry and GPU_Geometry classes are defined in
	// Geometry.h/Geometry.cpp They will work for this assignment, but for some of
	// the bonuses you may have to modify them.
	//initialize geometry objects
	CPU_Geometry controlPtsCPUGeom;
	CPU_Geometry selectedControlPtsCPUGeom;
	CPU_Geometry curveCPUGeom;
	CPU_Geometry surfaceOfRevolutionPtsCPUGeom;
	CPU_Geometry surfaceOfRevolutionTriCPUGeom;
	CPU_Geometry exSurfCntrlPtsCPUGeom;
	CPU_Geometry exSurfCntrlTriCPUGeom;
	CPU_Geometry exSurfPtsCPUGeom;
	CPU_Geometry exSurfTriCPUGeom;
	GPU_Geometry controlPtsGPUGeom;
	GPU_Geometry selectedControlPtsGPUGeom;
	GPU_Geometry controlLinesGPUGeom;
	GPU_Geometry curveGPUGeom;
	GPU_Geometry surfaceOfRevolutionPtsGPUGeom;
	GPU_Geometry surfaceOfRevolutionTriGPUGeom;
	GPU_Geometry exSurfCntrlPtsGPUGeom;
	GPU_Geometry exSurfCntrlTriGPUGeom;
	GPU_Geometry exSurfPtsGPUGeom;
	GPU_Geometry exSurfTriGPUGeom;

	vector<glm::vec3> control_points;//initialize the control points vector
	vector<glm::vec3> selected_control_points;//initialise the selected control points vector

	//initialize controller peripheral state objects
	const MouseState& mouseState = callback_controller->getMouseState();
	const KeyboardState& keyboardState = callback_controller->getKeyboardState();

	glm::vec3* moveThisPoint;//pointer to a control point that is being clicked and dragged

	bool timeToUpdateReleaseLocation = false;

	// RENDER LOOP
	while (!window.shouldClose()) {
		glfwPollEvents();		

		//respond to mouse events
		if (callback_controller->mouseStateChanged()) {
			//update state
			if (mouseState.mouse_clicked_left && state.sceneNum == 1) {//if a user clicks to create a new points in 2D editor scene

				//get mouse coordinate pair
				//std::cout << "You clicked at: (" << mouseState.mouse_down_coord.x << ", " << mouseState.mouse_down_coord.y << ")" << std::endl;

				bool clickedExisting = false;
				float distance = -1.f;
				//check if user clicked on an existing point
				for (auto it = control_points.begin(); it != control_points.end(); it++) {//clicked on an existing point
					distance = glm::distance(glm::vec3(mouseState.mouse_down_coord, 0.f), *it);
					if (distance < 0.05 && distance > 0.0) {
						clickedExisting = true;
						callback_controller->setMovingPoint(true);
						moveThisPoint = &*it;
						//std::cout << "Clicked an existing point!" << std::endl;
					}
				}
				//place the coordinate pair into control points vector if it is not an existing point
				if (clickedExisting == false)
					control_points.push_back(glm::vec3(mouseState.mouse_down_coord, 0.f));
				state.pointChanged = true;
			}

			else if (mouseState.mouse_clicked_right && state.sceneNum == 1) {//if a user right clicks on a point then it is selected
				selectPointOnClick(mouseState, control_points, selected_control_points, state);
			}

			else if (mouseState.mouse_clicked_left && (state.sceneNum == 2 || state.sceneNum == 3 || state.sceneNum == 4)) {//store location of mouse click
				//get mouse coordinate pair where the left button is clicked
				camera_down_coord = mouseState.mouse_down_coord;
				callback_controller->setMovingCamera(true);
			}

			callback_controller->mouseStateHandled();//turn off once the state change is addressed because we don't want to update every frame, just when changes occur
		}

		else if (state.movingPoint == true && state.sceneNum == 1) {//drag and drop of existing point is performed in this loop
			moveExistingPoint(control_points, callback_controller, moveThisPoint, controlPtsCPUGeom, controlPtsGPUGeom);
		}

		else if (state.movingCamera == true && (state.sceneNum == 2 || state.sceneNum == 3 || state.sceneNum == 4)) {//rotate camera direction with mouse events
			MouseState tempMouseState = callback_controller->getMouseState();

			if (tempMouseState.mouse_held_left == true) {//update interface while left mouse button is pressed down
				camera_down_coord = mouseState.mouse_down_coord;
			}
			else if (tempMouseState.mouse_release_left == true) {//update the final location when the left mouse button is released
				callback_controller->setMovingCamera(false);
				camera_up_coord = mouseState.mouse_up_coord;

				//update camera direction
				cameraXoffset = -(camera_up_coord.x - camera_down_coord.x);
				cameraYoffset = (camera_down_coord.y - camera_up_coord.y);

				cameraXoffset = cameraXoffset * cameraSensitivity;
				cameraYoffset = cameraYoffset * cameraSensitivity;

				yaw = yaw + cameraXoffset;
				pitch = pitch + cameraYoffset;

				//to prevent bending the neck backwards and looking through our own body we limit the camera orientation similar to common FPS games
				if (pitch > 89.0f)
					pitch = 89.0f;
				if (pitch < -89.0f)
					pitch = -89.0f;

				glm::vec3 direction;
				direction.x = cos(glm::radians(yaw)) * cos(glm::radians(pitch));
				direction.y = sin(glm::radians(pitch));
				direction.z = sin(glm::radians(yaw)) * cos(glm::radians(pitch));
				cameraFront = glm::normalize(direction);
			}
		}

		//respond to keyboard events
		else if (callback_controller->keyboardStateChanged()) {
			if (keyboardState.delete_pressed && state.sceneNum == 1) {//on press of the delete key delete the currently selected points
				deleteSelectedPts(control_points, selected_control_points, state);
			}
			else if (keyboardState.r_pressed && state.sceneNum == 1) {//on press of the R key reset all objects that are used by the program
				resetPts(control_points, selected_control_points, controlPtsCPUGeom, selectedControlPtsCPUGeom, curveCPUGeom, surfaceOfRevolutionPtsCPUGeom);
			}
			else if (keyboardState.t_pressed && (state.sceneNum == 3 || state.sceneNum == 4)) {//on press of the T key toggle the wireframe and solid draw mode in scene 3
				toggleWireframeSolidView(state, callback_controller);
			}
			callback_controller->keyboardStatehandled();
		}
		
		//check if user changed scene parameters and update the cpuGeom objects accordingly
		if (callback_controller->getState() != state) {

			//update the state object inside main to match the state private variable of the callback
			state = callback_controller->getState();

			//clear the cpuGeom object vertices and their colors
			controlPtsCPUGeom.verts.clear();
			controlPtsCPUGeom.cols.clear();
			selectedControlPtsCPUGeom.verts.clear();
			selectedControlPtsCPUGeom.cols.clear();
			curveCPUGeom.verts.clear();
			curveCPUGeom.cols.clear();
			surfaceOfRevolutionPtsCPUGeom.verts.clear();
			surfaceOfRevolutionPtsCPUGeom.cols.clear();
			surfaceOfRevolutionTriCPUGeom.verts.clear();
			surfaceOfRevolutionTriCPUGeom.cols.clear();
			exSurfCntrlPtsCPUGeom.verts.clear();
			exSurfCntrlPtsCPUGeom.cols.clear();
			exSurfCntrlTriCPUGeom.verts.clear();
			exSurfCntrlTriCPUGeom.cols.clear();
			exSurfPtsCPUGeom.verts.clear();
			exSurfPtsCPUGeom.cols.clear();
			exSurfTriCPUGeom.verts.clear();
			exSurfTriCPUGeom.cols.clear();

			//Generate cpuGeom object that reflects the user selected scene
			if ((state.sceneNum == 1 || state.sceneNum == 2) && state.subSceneNum == 1 && control_points.size()>0) {//2D Bezier Edit Scene and 3D Bezier Viewer Scene
				setupBezierEditScene(state, control_points, controlPtsCPUGeom, controlPtsGPUGeom,
					controlLinesGPUGeom, curveCPUGeom, curveGPUGeom, selected_control_points,
					selectedControlPtsCPUGeom, selectedControlPtsGPUGeom);
			}

			else if ((state.sceneNum == 1 || state.sceneNum == 2) && state.subSceneNum == 2) {//2D B-Spline Edit scene and 3D B-Spline Viewer scene
				//pass the control points to generate the B-Spline curve
				if(control_points.size()>1)
				setupBSplineEditScene(state, control_points, controlPtsCPUGeom, controlPtsGPUGeom,
					controlLinesGPUGeom, curveCPUGeom, curveGPUGeom, selected_control_points,
					selectedControlPtsCPUGeom, selectedControlPtsGPUGeom);
				else if (control_points.size() == 1)
					setupBezierEditScene(state, control_points, controlPtsCPUGeom, controlPtsGPUGeom,
						controlLinesGPUGeom, curveCPUGeom, curveGPUGeom, selected_control_points,
						selectedControlPtsCPUGeom, selectedControlPtsGPUGeom);
			}

			else if (state.sceneNum == 3 && (state.subSceneNum == 1 || state.subSceneNum == 2)) {//3D surface of rotation viwer scene
				//set up the surface of revolution scene
				setupSurfaceOfRevolutionScene(control_points, controlPtsCPUGeom, controlPtsGPUGeom,
											  controlLinesGPUGeom, state, curveCPUGeom, curveGPUGeom,
											  selected_control_points, selectedControlPtsCPUGeom,
											  selectedControlPtsGPUGeom, surfaceOfRevolutionTriCPUGeom,
											  surfaceOfRevolutionTriGPUGeom);
			}
			else if (state.sceneNum == 4 && state.subSceneNum == 1 ) {//3D scene for a first Example of the Tensor Product Sufrace scene
				setUpTensorProductSurfaceEx1Scene(exSurfCntrlPtsCPUGeom, exSurfCntrlTriCPUGeom,	exSurfTriCPUGeom,
												  exSurfPtsCPUGeom, exSurfCntrlPtsGPUGeom, exSurfCntrlTriGPUGeom,
												  exSurfPtsGPUGeom, exSurfTriGPUGeom);
			}
			else if (state.sceneNum == 4 && state.subSceneNum == 2) {//3D scene for a second Example of the Tensor Product Sufrace scene
				setUpTensorProductSurfaceEx2Scene(exSurfCntrlPtsCPUGeom, exSurfCntrlTriCPUGeom, exSurfTriCPUGeom,
												  exSurfPtsCPUGeom, exSurfCntrlPtsGPUGeom, exSurfCntrlTriGPUGeom,
												  exSurfPtsGPUGeom, exSurfTriGPUGeom);
			}
		}

		glEnable(GL_LINE_SMOOTH);
		glEnable(GL_FRAMEBUFFER_SRGB);
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		shader.use();

		//Setup matricies for 3D movement of the camera and link them to our vertex shader.
		GLint locationOfModelMat = glGetUniformLocation(shader.getID(), "model");
		GLint locationOfViewMat = glGetUniformLocation(shader.getID(), "view");
		GLint locationOfProjectionMat = glGetUniformLocation(shader.getID(), "projection");
		glUniformMatrix4fv(locationOfModelMat, 1, false, glm::value_ptr(model));
		glUniformMatrix4fv(locationOfViewMat, 1, false, glm::value_ptr(view));
		glUniformMatrix4fv(locationOfProjectionMat, 1, false, glm::value_ptr(proj));

		//Select the draw methods to accomodate the current scene
		if (state.sceneNum == 1 && (state.subSceneNum == 1 || state.subSceneNum == 2)) {
			//update matricies used for 3D movement of the camera leave as Identity since no movement takes place in this scene
			model = glm::mat4(1.0f);
			view = glm::mat4(1.0f);
			proj = glm::mat4(1.0f);			

			//draw the lines between control points
			controlLinesGPUGeom.bind();
			glDrawArrays(GL_LINE_STRIP, 0, GLsizei(controlPtsCPUGeom.verts.size()));

			//draw the control points
			controlPtsGPUGeom.bind();
			glDrawArrays(GL_POINTS, 0, GLsizei(controlPtsCPUGeom.verts.size()));

			//draw the Bezier curve
			curveGPUGeom.bind();
			glDrawArrays(GL_LINE_STRIP, 0, GLsizei(curveCPUGeom.verts.size()));

			//draw the selected points
			selectedControlPtsGPUGeom.bind();
			glDrawArrays(GL_POINTS, 0, GLsizei(selectedControlPtsCPUGeom.verts.size()));
		}
		else if (state.sceneNum == 2 && (state.subSceneNum == 1 || state.subSceneNum == 2)) {			
			//update matricies used for 3D movement of the camera
			model = glm::mat4(1.0f);
			view = glm::lookAt(cameraPos, cameraPos + cameraFront, cameraUp);
			proj = glm::perspective(glm::radians(55.0f), (float)800.0 / (float)800.0, 0.1f, 100.0f);

			//draw the lines between control points
			controlLinesGPUGeom.bind();
			glDrawArrays(GL_LINE_STRIP, 0, GLsizei(controlPtsCPUGeom.verts.size()));

			//draw the control points
			controlPtsGPUGeom.bind();
			glDrawArrays(GL_POINTS, 0, GLsizei(controlPtsCPUGeom.verts.size()));

			//draw the Bezier curve
			curveGPUGeom.bind();
			glDrawArrays(GL_LINE_STRIP, 0, GLsizei(curveCPUGeom.verts.size()));

			//draw the selected points
			selectedControlPtsGPUGeom.bind();
			glDrawArrays(GL_POINTS, 0, GLsizei(selectedControlPtsCPUGeom.verts.size()));
		}
		else if (state.sceneNum == 3 && (state.subSceneNum == 1 || state.subSceneNum == 2)) {
			//update matricies used for 3D movement of the camera
			model = glm::mat4(1.0f);
			view = glm::lookAt(cameraPos, cameraPos + cameraFront, cameraUp);
			proj = glm::perspective(glm::radians(55.0f), (float)800.0 / (float)800.0, 0.1f, 100.0f);

			//draw the quads with triangles
			surfaceOfRevolutionTriGPUGeom.bind();
			glDrawArrays(GL_TRIANGLES, 0, (GLsizei)surfaceOfRevolutionTriCPUGeom.verts.size());
		}
		else if (state.sceneNum == 4 && (state.subSceneNum == 1|| state.subSceneNum == 2)) {
			//update matricies used for 3D movement of the camera
			model = glm::mat4(1.0f);
			view = glm::lookAt(cameraPos, cameraPos + cameraFront, cameraUp);
			proj = glm::perspective(glm::radians(55.0f), (float)800.0 / (float)800.0, 0.1f, 100.0f);

			//draw the control points
			exSurfCntrlPtsGPUGeom.bind();
			glDrawArrays(GL_POINTS, 0, GLsizei(exSurfCntrlPtsCPUGeom.verts.size()));

			//draw the control surface with triangles
			exSurfCntrlTriGPUGeom.bind();
			glDrawArrays(GL_TRIANGLES, 0, GLsizei(exSurfCntrlTriCPUGeom.verts.size()));

			//draw the tensor product sufrace with triangles
			exSurfTriGPUGeom.bind();
			glDrawArrays(GL_TRIANGLES, 0, GLsizei(exSurfTriCPUGeom.verts.size()));
		}

		glDisable(GL_FRAMEBUFFER_SRGB); // disable sRGB for things like imgui

		//Display scene information on the screen
		// Starting the new ImGui frame
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();
		// Putting the text-containing window in the top-left of the screen.
		ImGui::SetNextWindowPos(ImVec2(5, 5));
		// Setting flags
		ImGuiWindowFlags textWindowFlags =
			ImGuiWindowFlags_NoMove |				// text "window" should not move
			ImGuiWindowFlags_NoResize |				// should not resize
			ImGuiWindowFlags_NoCollapse |			// should not collapse
			ImGuiWindowFlags_NoSavedSettings |		// don't want saved settings mucking things up
			ImGuiWindowFlags_AlwaysAutoResize |	    // window should auto-resize to fit the text
			//ImGuiWindowFlags_NoBackground |		// window should be transparent; only the text should be visible
			ImGuiWindowFlags_NoDecoration |			// no decoration; only the text should be visible
			ImGuiWindowFlags_NoTitleBar;			// no title; only the text should be visible

		// Begin a new window with these flags. (bool *)0 is the "default" value for its argument.
		ImGui::Begin("scoreText", (bool*)0, textWindowFlags);

		// Scale up text a little, and set its value
		ImGui::SetWindowFontScale(1.5f);
		if (state.sceneNum == 1 && state.subSceneNum == 1) {
			ImGui::Text("Scene Name: 2D Bezier Curve Editor");
			ImGui::Text("Scene Number: 1.1");
			//ImGui::Text("Sub-scene Number: 1");
		}
		else if (state.sceneNum == 1 && state.subSceneNum == 2) {
			ImGui::Text("Scene Name: 2D B-Spline Curve Editor");
			ImGui::Text("Scene Number: 1.2");
			//ImGui::Text("Sub-scene Number: 2");
		}
		else if (state.sceneNum == 2 && state.subSceneNum == 1) {
			ImGui::Text("Scene Name: 3D Bezier Curve Viewer");
			ImGui::Text("Scene Number: 2.1");
			//ImGui::Text("Sub-scene Number: 1");
		}
		else if (state.sceneNum == 2 && state.subSceneNum == 2) {
			ImGui::Text("Scene Name: 3D B-Spline Curve Viewer");
			ImGui::Text("Scene Number: 2.2");
			//ImGui::Text("Sub-scene Number: 2");
		}
		else if (state.sceneNum == 3) {
			ImGui::Text("Scene Name: 3D B-Spline Surface of Revolution Viewer");
			ImGui::Text("Scene Number: 3");
		}
		else if (state.sceneNum == 4 && state.subSceneNum == 1) {
			ImGui::Text("Scene Name: 3D Bezier Tensor Product Surface Example");
			ImGui::Text("Scene Number: 4.1");
			//ImGui::Text("Sub-scene Number: 1");
		}
		else if (state.sceneNum == 4 && state.subSceneNum == 2) {
			ImGui::Text("Scene Name: 3D B-Spline Tensor Product Surface Example");
			ImGui::Text("Scene Number: 4.2");
			//ImGui::Text("Sub-scene Number: 2");
		}
		// End the window.
		ImGui::End();
		ImGui::Render();	// Render the ImGui window
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData()); // Some middleware thing

		window.swapBuffers();//make new things appear on the screen
	}
	// ImGui cleanup
	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();

	glfwTerminate();
	return 0;
}
