#include <algorithm>
#include <format>
#include <iostream>
#include <vector>
#include <list>

#include <imgui.h>
#include <misc/cpp/imgui_stdlib.h>
#include <backends/imgui_impl_glfw.h>
#include <backends/imgui_impl_opengl2.h>
#include <GLFW/glfw3.h>

#define LOG_ERROR(text, ...)	std::cout << std::format(text, __VA_ARGS__)

import MassHelper;
import PeriodicTable;

import UtlWinConsole;

using std::vector;
using std::string;
using std::list;

// Example 1 SAMPLER
//vector<MassPeak_t> g_rgflMassData = { 86, 113, 131, 141, 159, 175, 158, 262, 286, 290, 304, 387, 402, 417, 500, 514, 611, 597, 629, 645, 716 };
//double g_flMPlusOne = 402 * 2 - 1;
//string g_szInputMassData = "86\n113\n131\n141\n159\n175\n158\n262\n286\n290\n304\n387\n402\n417\n500\n514\n611\n597\n629\n645\n716";

// Example 2 CHEMISTK
//vector<MassPeak_t> g_rgflMassData = { 86, 110, 147, 207, 230, 248, 270, 298, 317, 335, 409, 399, 427, 448, 503, 510, 558, 579, 671, 708, 845 };
//double g_flMPlusOne = 503 * 2 - 1;
//string g_szInputMassData = "86\n110\n147\n207\n230\n248\n270\n298\n317\n335\n409\n399\n427\n448\n503\n510\n558\n579\n671\n708\n845";

// Example 3 DIERNIVAGK
//vector<MassPeak_t> g_rgflMassData = { 86.0, 129.0, 147.1, 201.1, 204.1, 229.1, 275.1, 358.1, 347.2, 379.2, 443.7, 484.5, 487.3, 500.3, 514.2, 557.5, 497.2, 601.3668, 628.3, 741.4, 757.4, 840.4, 886.5, 911.5, 999.6 };
//double g_flMPlusOne = 1114.6215;
//string g_szInputMassData = "86.0\n129.0\n147.1\n201.1\n204.1\n229.1\n275.1\n358.1\n347.2\n379.2\n443.7\n484.5\n487.3\n500.3\n514.2\n557.5\n497.2\n601.3668\n628.3\n741.4\n757.4\n840.4\n886.5\n911.5\n999.6";

// Example 4 IFSQVGK
//vector<MassPeak_t> g_rgflMassData = { 128.08, 147.11, 204.13, 233.16, 261.16, 303.20, 348.19, 380.72, 389.73, 431.26, 458.24, 476.25, 500.28, 518.29, 557.31, 575.32, 632.34, 647.35, 665.36 };
//double g_flMPlusOne = 778.4458;
//string g_szInputMassData = "128.08\n147.11\n204.13\n233.16\n261.16\n303.20\n348.19\n380.72\n389.73\n431.26\n458.24\n476.25\n500.28\n518.29\n557.31\n575.32\n632.34\n647.35\n665.36";

// Example 5 HGIWNYK
//vector<MassPeak_t> g_rgflMassData = { 308.2, 310.2, 424.2, 459.2, 494.3, 608.3, 610.3, 723.4, 771.4, 780.4 };
//double g_flMPlusOne = 459.2 * 2 - 1;
//string g_szInputMassData = "308.2\n310.2\n424.2\n459.2\n494.3\n608.3\n610.3\n723.4\n771.4\n780.4";

// Example 6 NTDGTQIIGYmTVNSR
vector<MassPeak_t> g_rgflMassData = { 548.2, 558.2, 576.2, 617.3, 659.3, 723.3, 730.4, 843.5, 886.3, 893.4, 900.5, 943.3, 1056.4, 1063.6, 1146.6, 1192.6, 1210.6, 1169.5, 1297.6, 1311.7, 1380.7, 1398.7, 1410.7, 1455.7, 1524.7, 1552.7, 1570.8, 1611.7, 1671.8 };
double g_flMPlusOne = 1785.8436;
string g_szInputMassData = "548.2\n558.2\n576.2\n617.3\n659.3\n723.3\n730.4\n843.5\n886.3\n893.4\n900.5\n943.3\n1056.4\n1063.6\n1146.6\n1192.6\n1210.6\n1169.5\n1297.6\n1311.7\n1380.7\n1398.7\n1410.7\n1455.7\n1524.7\n1552.7\n1570.8\n1611.7\n1671.8";

double g_flMPlusTwo = (g_flMPlusOne + amu::Hydrogen) / 2.0;

list<AlternativeReality_t> g_rgExplanations;
string g_szSelectedWorldline;

void DrawInputWindow(void) noexcept
{
	bool bChanged = false;

	ImGui::Begin("Input");

	if (ImGui::InputDouble("[M+1] ion mass", &g_flMPlusOne))
	{
		bChanged = true;
		g_flMPlusTwo = (g_flMPlusOne + amu::Hydrogen) / 2.0;
	}

	if (ImGui::InputDouble("[M+2] ion mass", &g_flMPlusTwo))
	{
		bChanged = true;
		g_flMPlusOne = g_flMPlusTwo * 2 - amu::Hydrogen;
	}

	if (ImGui::InputTextMultiline("Peaks", &g_szInputMassData))
	{
		bChanged = true;
		g_rgflMassData.clear();

		auto lastPos = g_szInputMassData.find_first_not_of(",\n \t", 0);
		auto pos = g_szInputMassData.find_first_of(",\n \t", lastPos);

		while (string::npos != pos || string::npos != lastPos)
		{
			g_rgflMassData.emplace_back(std::stod(g_szInputMassData.substr(lastPos, pos - lastPos)));
			lastPos = g_szInputMassData.find_first_not_of("\n \t", pos);
			pos = g_szInputMassData.find_first_of("\n \t", lastPos);
		}

		std::sort(g_rgflMassData.begin(), g_rgflMassData.end(), std::less<MassPeak_t>());
	}

	if (ImGui::Button("Deduce") && g_flMPlusOne > 0)
	{
		//std::cout << "\n\nNew analysis begin:\n";
		//IdentifyBorderIons(g_rgflMassData, M_plus_1);
		//RecursiveIdentify(g_rgflMassData, IonType::b, -2);
		//RecursiveIdentify(g_rgflMassData, IonType::y, -2);
		//ParseSpectrum<IonType::y>(g_rgflMassData, M_plus_1);
		//ParseSpectrum<IonType::b>(g_rgflMassData, M_plus_1);

		clear_console();
		g_rgExplanations = Solve(g_rgflMassData, g_flMPlusOne);
	}

	ImGui::NewLine();
	ImGui::NewLine();

	static MassPeak_t SelectedPeak{};

	if (ImGui::BeginTable("Result", 4, ImGuiTableFlags_Resizable))
	{
		// Showing the '0' option.
		ImGui::TableNextRow();
		ImGui::TableSetColumnIndex(0);
		if (ImGui::RadioButton("Naught", SelectedPeak == 0.0))
			SelectedPeak = MassPeak_t{};
		ImGui::TableSetColumnIndex(1);
		ImGui::TextUnformatted(std::to_string(-SelectedPeak.m_Value).c_str());
		ImGui::TableSetColumnIndex(3);
		ImGui::TextUnformatted(TestNumber3(SelectedPeak.m_Value).c_str());

		// Showing actually peak.
		for (const auto& Peak : g_rgflMassData)
		{
			ImGui::TableNextRow();

			ImGui::TableSetColumnIndex(0);
			if (ImGui::RadioButton(std::to_string(Peak.m_Value).c_str(), SelectedPeak == Peak))
				SelectedPeak = Peak;
			ImGui::TableSetColumnIndex(1);
			ImGui::TextUnformatted(std::to_string(Peak - SelectedPeak).c_str());
			ImGui::TableSetColumnIndex(2);
			ImGui::TextUnformatted(Peak.ToString().c_str());
			ImGui::TableSetColumnIndex(3);
			ImGui::TextUnformatted(TestNumber3(std::abs(SelectedPeak - Peak)).c_str());
		}

		// Showing the [M+1] option.
		ImGui::TableNextRow();
		ImGui::TableSetColumnIndex(0);
		if (ImGui::RadioButton("[M+1]", SelectedPeak == g_flMPlusOne))
			SelectedPeak = g_flMPlusOne;
		ImGui::TableSetColumnIndex(1);
		ImGui::TextUnformatted(std::to_string(g_flMPlusOne - SelectedPeak).c_str());
		ImGui::TableSetColumnIndex(3);
		ImGui::TextUnformatted(TestNumber3(std::abs(SelectedPeak - g_flMPlusOne)).c_str());

		ImGui::EndTable();
	}

	ImGui::End();
}

void DrawAnalyzeWindow(void) noexcept
{
	ImGui::Begin("Analyze");

	ImGui::Bullet(); ImGui::SameLine(); ImGui::TextUnformatted(std::format("[M+H]: {}", g_flMPlusOne).c_str());
	ImGui::Bullet(); ImGui::SameLine(); ImGui::TextUnformatted(std::format("[M+H] - H2O: {}", g_flMPlusOne - amu::Hydrogen * 2 - amu::Oxygen).c_str());
	ImGui::Bullet(); ImGui::SameLine(); ImGui::TextUnformatted(std::format("[M+H] - NH3: {}", g_flMPlusOne - amu::Hydrogen * 3 - amu::Nitrogen).c_str());

	for (const auto& Worldline : g_rgExplanations)
	{
		auto szSeq = Conclude(Worldline.m_Solution);
		if (ImGui::RadioButton(szSeq.c_str(), g_szSelectedWorldline == szSeq))
		{
			g_szSelectedWorldline = szSeq;

			ResetPeaks(g_rgflMassData);
			MarkPeaks(Worldline, g_rgflMassData);
		}
	}

	ImGui::End();
}

int main(int, char**) noexcept
{
	// 101 basically setup for ANY C++ project.
	std::ios_base::sync_with_stdio(false);

	// Setup window
	glfwSetErrorCallback([](int error, const char* description) { LOG_ERROR("Error {}: {}\n", error, description); });

	if (!glfwInit())
		return EXIT_FAILURE;

	GLFWwindow* window = glfwCreateWindow(1280, 720, "Mass Helper", NULL, NULL);
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

		// Input window.
		DrawInputWindow();
		DrawAnalyzeWindow();

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
