/**
 * @file
 *
 * Implementation of useful test routines.
 *
 * @ingroup relic
 */
#include "test.h"

/*
 * Turn on colors for pretty printing.
 */
#define COLOR

#ifdef _MSC_VER
#undef DOUBLE
#include <Windows.h>
/**
 * Stores the default color used in the console.
 */
static int default_color;
#endif

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

/**
 * Color of the string printed when the test fails (red).
 */
#ifdef _MSC_VER
#define FAIL_COLOR      12
#else
#define FAIL_COLOR		31
#endif

/**
 * Color of the string printed when the test passes (green).
 */
#ifdef _MSC_VER
#define PASS_COLOR      10
#else
#define PASS_COLOR		32
#endif

/**
 * Command to set terminal colors.
 */
#define CMD_SET			27

/**
 * Command to reset terminal colors.
 */
#define CMD_RESET		0

/**
 * Print with bright attribute.
 */
#define CMD_ATTR		1

/**
 * Copies default color to global variable.
 */
static void cache_default_color(void) {
#ifdef _MSC_VER
	CONSOLE_SCREEN_BUFFER_INFO csbi;
	HANDLE m_hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
	GetConsoleScreenBufferInfo(m_hConsole, &csbi);
	default_color = (csbi.wAttributes & 255);
#endif
}

/**
 * Changes font to test failure mode.
 */
static void fail_font(void) {
    cache_default_color();
#ifdef COLOR
#ifdef _MSC_VER
	HANDLE m_hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
	SetConsoleTextAttribute(m_hConsole, FAIL_COLOR);
#else
	printf("%c[%d;%dm", CMD_SET, CMD_ATTR, FAIL_COLOR);
#endif
#endif
}

/**
 * Changes font to test pass mode.
 */
static void pass_font(void) {
    cache_default_color();
#ifdef COLOR
#ifdef _MSC_VER
	HANDLE m_hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
	SetConsoleTextAttribute(m_hConsole, PASS_COLOR);
#else
	printf("%c[%d;%dm", CMD_SET, CMD_ATTR, PASS_COLOR);
#endif
#endif
}

/**
 * Resets font to default.
 */
static void reset_font(void) {
#ifdef COLOR
#ifdef _MSC_VER
	HANDLE m_hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
	SetConsoleTextAttribute(m_hConsole, default_color);
#else
	printf("%c[%dm", CMD_SET, CMD_RESET);
#endif
#endif
}

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

void test_fail(void) {
	fail_font();
	printf("[FAIL]\n");
	reset_font();
}

void test_pass(void) {
	pass_font();
	printf("[PASS]\n");
	reset_font();
}
