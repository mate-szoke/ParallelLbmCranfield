/**
 * Macros to manipulate boundary conditions bitmask
 * @file BcMacros.h
 * @author Adam Koleszar (adam.koleszar@gmail.com)
 * @details The bitmask looks the following way
 * @verbatim
   | 30 | 28 | 26 | 24 | 22 | 20 | 18 | 16 | 14 | 12 | 10 |  8 |  6 |  4 |  2 |  0 |
   |-------------------|---------|----|----|----|----|----|----|----|----|----|----|
   | boundary ID (8bit)|  empty  | FC | M0 | SE | SW | NW | NE |  S |  W |  N |  E |@endverbatim
 */
#ifndef BC_MACROS_H
#define BC_MACROS_H

#define BC_WALL_E   0x00001 ///< BC: east  (1) wall
#define BC_INLT_E   0x00002 ///< BC: east  (1) inlet
#define BC_OUTL_E   0x00003 ///< BC: east  (1) outlet
#define BC_WALL_N   0x00004 ///< BC: north (2) wall
#define BC_INLT_N   0x00008 ///< BC: north (2) inlet
#define BC_OUTL_N   0x0000C ///< BC: north (2) outlet
#define BC_WALL_W   0x00010 ///< BC: west  (3) wall
#define BC_INLT_W   0x00020 ///< BC: west  (3) inlet
#define BC_OUTL_W   0x00030 ///< BC: west  (3) outlet
#define BC_WALL_S   0x00040 ///< BC: south (4) wall
#define BC_INLT_S   0x00080 ///< BC: south (4) inlet
#define BC_OUTL_S   0x000C0 ///< BC: south (4) outlet

#define BC_WALL_NE  0x00100 ///< BC: north-east (5) wall
#define BC_INLT_NE  0x00200 ///< BC: north-east (5) inlet
#define BC_OUTL_NE  0x00300 ///< BC: north-east (5) outlet
#define BC_WALL_NW  0x00400 ///< BC: north-west (6) wall
#define BC_INLT_NW  0x00800 ///< BC: north-west (6) inlet
#define BC_OUTL_NW  0x00C00 ///< BC: north-west (6) outlet
#define BC_WALL_SW  0x01000 ///< BC: south-west (7) wall
#define BC_INLT_SW  0x02000 ///< BC: south-west (7) inlet
#define BC_OUTL_SW  0x03000 ///< BC: south-west (7) outlet
#define BC_WALL_SE  0x04000 ///< BC: south-east (8) wall
#define BC_INLT_SE  0x08000 ///< BC: south-east (8) inlet
#define BC_OUTL_SE  0x0C000 ///< BC: south-east (8) outlet

#define BC_WALL_B   0x10000 ///< BC: main (0) wall
#define BC_INLT_B   0x20000 ///< BC: main (0) inlet
#define BC_OUTL_B   0x30000 ///< BC: main (0) oulet
#define BC_CORNER   0x40000 ///< BC: corner
#define BC_FLUID    0x80000 ///< BC: fluid

#define BC_NONE     0x0 ///< BC: type none    (00)
#define BC_WALL     0x1 ///< BC: type wall    (01)
#define BC_INLT     0x2 ///< BC: type inlet   (10)
#define BC_OUTL     0x3 ///< BC: type oulet   (11)
#define BC_ALL      0x3 ///< BC: type all/any (11)

#define BC_E(i)  ((i)<< 0)    ///< BC: set type to east       @param i BC type
#define BC_N(i)  ((i)<< 2)    ///< BC: set type to north      @param i BC type
#define BC_W(i)  ((i)<< 4)    ///< BC: set type to west       @param i BC type
#define BC_S(i)  ((i)<< 6)    ///< BC: set type to south      @param i BC type
#define BC_NE(i) ((i)<< 8)    ///< BC: set type to north-east @param i BC type
#define BC_NW(i) ((i)<<10)    ///< BC: set type to north-west @param i BC type
#define BC_SW(i) ((i)<<12)    ///< BC: set type to south-west @param i BC type
#define BC_SE(i) ((i)<<14)    ///< BC: set type to south-east @param i BC type
#define BC_B(i)  ((i)<<16)    ///< BC: set type to main       @param i BC type
#define BC_C(i)  ((i)<<18)    ///< BC: set type to corner     @param i true/false
#define BC_F(i)  ((i)<<19)    ///< BC: set type to fluid      @param i true/false

/**
 * @brief BC: set type to specified direction
 *
 * @param type BC type
 * @param dir direction
 * @return bitmask
 */
#define BC_MASK(type, dir) ( (type) << (((dir)-1) * 2) )

#define BND_ID_ALL 0xFF000000 ///< BC: boundary id all
#define BOUND_ID(i) ((i)<<24) ///< BC: set boundary id for lift/drag @param i boundary ID

#endif