/**
 * Check CUDA API calls
 * @file Check.h
 * @author Adam Koleszar (adam.koleszar@gmail.com)
 */
#ifndef CHECK_H
#define CHECK_H

/**
 * @brief Check return value of CUDA API calls
 *
 * @param err return value from API call
 * @param location source file where check was called
 * @param lineno line number where check was called
 * @warning exits program if error occurs without cleanup
 */
void check(cudaError err, const char *location, int lineno);

/**
 * @brief Check macro for easier access
 *
 * @param err return value from API call
 */
#define CHECK(err) check(err, __FILE__, __LINE__)

#endif
