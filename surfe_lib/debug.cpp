#include <debug.h>

void open_console_window()
{
	// Grab a console to show debug info ...
	AllocConsole();
	freopen("CONOUT$","w",stdout);

	HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
	_CONSOLE_SCREEN_BUFFER_INFO info;
	GetConsoleScreenBufferInfo(hConsole,&info);
	SMALL_RECT windowMaxSize = {0,0,info.dwMaximumWindowSize.X-1,info.dwMaximumWindowSize.Y-1};
	SetConsoleWindowInfo(hConsole,TRUE,&windowMaxSize);
	COORD coord;
	coord.X = 3000;
	coord.Y = 8000;
	SetConsoleScreenBufferSize(hConsole,coord);
	SetConsoleTitle("Generalized Interpolation with Radial Basis Functions");
}
