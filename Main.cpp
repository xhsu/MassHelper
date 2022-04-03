#include <iostream>
#include <format>
#include <vector>

#include <imgui.h>
#include <misc/cpp/imgui_stdlib.h>
#include <backends/imgui_impl_glfw.h>
#include <backends/imgui_impl_opengl2.h>
#include <GLFW/glfw3.h>

#define LOG_ERROR(text, ...)	std::cout << std::format(text, __VA_ARGS__)

import MassHelper;

using std::vector;
using std::string;

vector<MassPeak_t> g_rgflMassData = { 308.2, 310.2, 424.2, 459.2, 494.3, 608.3, 610.3, 723.4, 771.4, 780.4 };
double M_plus_1 = 459.2 * 2 - 1;
// 803
string g_szInputMassData = "308.2, 310.2, 424.2, 459.2, 494.3, 608.3, 610.3, 723.4, 771.4, 780.4";
// 86, 113, 131, 141, 159, 175, 158, 262, 286, 290, 304, 387, 402, 417, 500, 514, 611, 597, 629, 645, 716

int main(int, char**) noexcept
{
	// Setup window
	glfwSetErrorCallback([](int error, const char* description) { LOG_ERROR("Error {}: {}\n", error, description); });

	if (!glfwInit())
		return EXIT_FAILURE;

	GLFWwindow* window = glfwCreateWindow(1280, 720, "Dear ImGui GLFW+OpenGL2 example", NULL, NULL);
	if (!window)
		return EXIT_FAILURE;

	glfwMakeContextCurrent(window);
	glfwSwapInterval(1); // Enable vsync

	// Setup Dear ImGui context
	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGuiIO& io = ImGui::GetIO(); (void)io;

	// Setup Dear ImGui style
	ImGui::StyleColorsClassic();

	// Setup Platform/Renderer backends
	ImGui_ImplGlfw_InitForOpenGL(window, true);
	ImGui_ImplOpenGL2_Init();

	// Main loop
	while (!glfwWindowShouldClose(window))
	{
		bool bChanged = false;

		// Poll and handle events (inputs, window resize, etc.)
		// You can read the io.WantCaptureMouse, io.WantCaptureKeyboard flags to tell if dear imgui wants to use your inputs.
		// - When io.WantCaptureMouse is true, do not dispatch mouse input data to your main application, or clear/overwrite your copy of the mouse data.
		// - When io.WantCaptureKeyboard is true, do not dispatch keyboard input data to your main application, or clear/overwrite your copy of the keyboard data.
		// Generally you may always pass all inputs to dear imgui, and hide them from your application based on those two flags.
		glfwPollEvents();

		// Start the Dear ImGui frame
		ImGui_ImplOpenGL2_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		ImGui::ShowDemoWindow();

		// 2. Show a simple window that we create ourselves. We use a Begin/End pair to created a named window.
		ImGui::Begin("Mass Helper");

		if (ImGui::InputDouble("[M+1] ion mass", &M_plus_1))
			bChanged = true;

		if (ImGui::InputTextMultiline("Peaks", &g_szInputMassData))
		{
			bChanged = true;
			g_rgflMassData.clear();

			auto lastPos = g_szInputMassData.find_first_not_of("\n \t", 0);
			auto pos = g_szInputMassData.find_first_of("\n \t", lastPos);

			while (string::npos != pos || string::npos != lastPos)
			{
				g_rgflMassData.emplace_back(std::stod(g_szInputMassData.substr(lastPos, pos - lastPos)));
				lastPos = g_szInputMassData.find_first_not_of("\n \t", pos);
				pos = g_szInputMassData.find_first_of("\n \t", lastPos);
			}
		}

		if (bChanged && M_plus_1 > 0)
		{
			std::cout << "\n\nNew analysis begin:\n";
			IdentifyBorderIons(g_rgflMassData, M_plus_1);
			RecursiveIdentify(g_rgflMassData, IonType::b, -2);
			RecursiveIdentify(g_rgflMassData, IonType::y, -2);
			ParseSpectrum<IonType::y>(g_rgflMassData, M_plus_1);
			ParseSpectrum<IonType::b>(g_rgflMassData, M_plus_1);
		}

		ImGui::NewLine();
		ImGui::NewLine();

		if (ImGui::BeginTable("Result", 2))
		{
			for (const auto& Peak : g_rgflMassData)
			{
				ImGui::TableNextRow();

				ImGui::TableSetColumnIndex(0);
				ImGui::TextUnformatted(std::to_string(Peak.m_Value).c_str());
				ImGui::TableSetColumnIndex(1);
				ImGui::TextUnformatted(Peak.ToString().c_str());

			}

			ImGui::EndTable();
		}

		ImGui::End();

		// Rendering
		ImGui::Render();
		int display_w, display_h;
		glfwGetFramebufferSize(window, &display_w, &display_h);
		glViewport(0, 0, display_w, display_h);
		glClearColor(0.45f, 0.55f, 0.60f, 1);
		glClear(GL_COLOR_BUFFER_BIT);

		// If you are using this code with non-legacy OpenGL header/contexts (which you should not, prefer using imgui_impl_opengl3.cpp!!),
		// you may need to backup/reset/restore other state, e.g. for current shader using the commented lines below.
		//GLint last_program;
		//glGetIntegerv(GL_CURRENT_PROGRAM, &last_program);
		//glUseProgram(0);
		ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
		//glUseProgram(last_program);

		glfwMakeContextCurrent(window);
		glfwSwapBuffers(window);
	}

	// Cleanup
	ImGui_ImplOpenGL2_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();

	glfwDestroyWindow(window);
	glfwTerminate();

	return EXIT_SUCCESS;
}
