/**
 * @defgroup tests Automated tests
 */

/**
 * @file
 *
 * Interface of useful routines for testing.
 *
 * @ingroup test
 */

#ifndef TEST_H
#define TEST_H

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

/*============================================================================*/
/* Macro definitions                                                          */
/*============================================================================*/

/**
 * Number of times a test is executed.
 */
#define TESTS 		100

/**
 * Runs a new benchmark once.
 *
 * @param[in] P				- the property description.
 */
#define TEST_ONCE(P)														\
	printf("Testing if " P "...%*c", (int)(64 - strlen(P)), ' ');			\

/**
 * Tests a sequence of commands to see if they respect some property.
 *
 * @param[in] P				- the property description.
 */
#define TEST_BEGIN(P)														\
	printf("Testing if " P "...%*c", (int)(64 - strlen(P)), ' ');			\
	for (int _i = 0; _i < TESTS; _i++)										\

/**
 * Asserts a condition.
 *
 * If the condition is not satisfied, a unconditional jump is made to the passed
 * label.
 *
 * @param[in] C				- the condition to assert.
 * @param[in] LABEL			- the label to jump if the condition is no satisfied.
 */
#define TEST_ASSERT(C, LABEL)												\
	if (!(C)) {																\
		test_fail();														\
		printf("(at ");														\
		printf(__FILE__);													\
		printf(":%d)\n", __LINE__);											\
		goto LABEL;															\
	}																		\

/**
 * Finalizes a test printing the test result.
 */
#define TEST_END															\
	test_pass()																\

/*============================================================================*/
/* Function prototypes                                                        */
/*============================================================================*/

/**
 * Prints a string indicating that the test failed.
 */
void test_fail(void);

/**
 * Prints a string indicating that the test passed.
 */
void test_pass(void);

#endif /* !RLC_TEST_H */
