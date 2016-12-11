/**
 * File to select floating point precision.
 *
 * Usage: make FLOAT_TYPE=USE_DOUBLE \<target\>
 *
 * Default: USE_SINGLE
 * @file FloatType.h
 * @author Adam Koleszar (adam.koleszar@gmail.com)
 */
#ifndef USE_DOUBLE
#ifndef FLOAT_TYPE
#define FLOAT_TYPE float    ///< floating point type (float/double)
#define FLOAT_FORMAT "%f"   ///< floating point formatting string for scanf
#endif
#else
#ifndef FLOAT_TYPE
#define FLOAT_TYPE double
#define FLOAT_FORMAT "%lf"
#endif
#endif
