/**
 * Command line parameter handling
 * @file Arguments.cu
 * @author Adam Koleszar (adam.koleszar@gmail.com)
 */
#include <string.h>
#if defined(_WIN32) || defined(_WIN64)
#include <regex>
#else
#include <regex.h>
#endif

#include "Arguments.h"
#include "FilesReading.h"
#include "argtable2.h"

/**
 * @brief Get inlet profile from command line parameter
 *
 * @param input command line parameter
 * @param defaultValue default value
 *
 * @return chosen inlet profile
 */
InletProfile getInletProfile(const char *input, InletProfile defaultValue)
{
    InletProfile ipr = defaultValue;
    if (strcmp(input, "yes") == 0)
    {
        ipr = INLET;
    }
    else if (strcmp(input, "no") == 0)
    {
        ipr = NO_INLET;
    }
    else if (strcmp(input, "pulsatile") == 0)
    {
        ipr = PULSATILE_INLET;
    }
    else if (strlen(input))
    {
        fprintf(stderr, "Invalid inlet profile (%s), going with default\n", input);
    }
    return ipr;
}

/**
 * @brief Get outlet profile from command line parameter
 *
 * @param input command line parameter
 * @param defaultValue default value
 *
 * @return chosen outlet profile
 */
OutletProfile getOutletProfile(const char *input, OutletProfile defaultValue)
{
    OutletProfile opr = defaultValue;
    if (!strcmp(input, "yes"))
    {
        opr = OUTLET;
    }
    else if (!strcmp(input, "second"))
    {
        opr = OUTLET_SECOND;
    }
    else if (!strcmp(input, "first"))
    {
        opr = OUTLET_FIRST;
    }
    else if (strlen(input))
    {
        fprintf(stderr, "Invalid outlet profile (%s), going with default\n", input);
    }
    return opr;
}

/**
 * @brief Get collision model from command line parameter
 *
 * @param input command line parameter
 * @param defaultValue default value
 *
 * @return chosen collision model
 */
CollisionModel getCollisionModel(const char *input, CollisionModel defaultValue)
{
    CollisionModel coll = defaultValue;
    if (!strcmp(input, "BGKW"))
    {
        coll = BGKW;
    }
    else if (!strcmp(input, "TRT"))
    {
        coll = TRT;
    }
    else if (!strcmp(input, "MRT"))
    {
        coll = MRT;
    }
    else if (strlen(input))
    {
        fprintf(stderr, "Invalid collision model (%s), going with default\n", input);
    }
    return coll;
}

/**
 * @brief Get boundary type from command line parameter
 *
 * @param input command line parameter
 * @param defaultValue default value
 *
 * @return chosen boundary type
 */
BoundaryType getBoundaryType(int input, BoundaryType defaultValue)
{
    BoundaryType bt = defaultValue;
    if (input > 0)
    {
        bt = CURVED;
    }
    else
    {
        bt = STRAIGHT;
    }
    return bt;
}

/**
 * @brief Get output format from command line parameter
 *
 * @param input command line parameter
 * @param defaultValue default value
 *
 * @return chosen output format
 */
OutputFormat getOutputFormat(const char *input, OutputFormat defaultValue)
{
    OutputFormat of = defaultValue;
    if (!strcmp(input, "paraview"))
    {
        of = PARAVIEW;
    }
    else if (!strcmp(input, "tecplot"))
    {
        of = TECPLOT;
    }
    else if (strlen(input))
    {
        fprintf(stderr, "Invalid output format (%s), going with default\n", input);
    }
    return of;
}

ArgResult handleArguments(int argc, char **argv, InputFilenames *inFn, Arguments *args)
{
    struct arg_lit  *helpArg  = arg_lit0 ("h", "help", "Print help options");
    struct arg_file *initArg  = arg_file0("f", "initfile", "<file>", "Initialisation from file (default: SetUpData.ini)");
    struct arg_lit  *testArg  = arg_lit0 ("t", "test", "Run unit tests");
    struct arg_file *nodeArg  = arg_file0("n", "node", "<file>", "Node file (default: Mesh/D2node.dat)");
    struct arg_file *bcArg = arg_file0("b", "bc", "<file>", "Boundary conditions file (default: Mesh/BCconnectors.dat)");
    struct arg_file *outArg   = arg_file0("o", "output", "<file>", "Output directory (default: Results)");

    struct arg_dbl *uArg      = arg_dbl0("u", "uavg", "<u>", "Mean U (x velocity)");
    struct arg_dbl *vArg      = arg_dbl0("v", "vavg", "<v>", "Mean V (y velocity)");
    struct arg_dbl *rhoArg    = arg_dbl0("r", "rho", "<rho>", "Density");
    struct arg_dbl *visArg    = arg_dbl0("s", "viscosity", "<nu>", "Viscosity");
#if defined(_WIN32) || defined(_WIN64)
    struct arg_str *inletArg  = arg_str0(NULL, "inlet", "[yes|no|pulsatile]", "Inlet profile (default: no)");
    struct arg_str *collArg   = arg_str0("c", "collision", "[BGKW|TRT|MRT]", "Collision model (default: BGKW)");
    struct arg_str *outlArg   = arg_str0("l", "outlet", "[yes|second|first]", "Outlet profile (default: second)");
    struct arg_str *formatArg = arg_str0(NULL, "format", "[paraview|tecplot]", "Output format (default: paraview)");
#else
    struct arg_rex *inletArg  = arg_rex0(NULL, "inlet", "[yes|no|pulsatile]", NULL, REG_ICASE, "Inlet profile (default: no)");
    struct arg_rex *collArg   = arg_rex0("c", "collision", "[BGKW|TRT|MRT]", NULL, REG_ICASE, "Collision model (default: BGKW)");
    struct arg_rex *outlArg   = arg_rex0("l", "outlet", "[yes|second|first]", NULL, REG_ICASE, "Outlet profile (default: second)");
    struct arg_rex *formatArg = arg_rex0(NULL, "format", "[paraview|tecplot]", NULL, REG_ICASE, "Output format (default: paraview)");
#endif
    struct arg_lit *curvedArg = arg_lit0(NULL, "curved", "Curved boundaries");
    struct arg_int *iterArg   = arg_int0("i", "iter", "<N>", "Number of iterations (default: 1000)");
    struct arg_int *autoeArg  = arg_int0(NULL, "every", "<N>", "Autosave after every <N> iterations (default: 0)");
    struct arg_int *autoaArg  = arg_int0(NULL, "after", "<N>", "Start autosaving after the <N>th iteration (default: 1000)");
    struct arg_int *drlftArg  = arg_int0("d", "draglift", "<id>", "Calculate drag/lift on <id> boundary (default: 0)");
    struct arg_end *endArg1   = arg_end (10);

    void *normalArgTable[] = {helpArg, initArg, testArg, nodeArg, bcArg, outArg, uArg, vArg, rhoArg, visArg,
                              inletArg, collArg, curvedArg, outlArg, iterArg, autoeArg, autoaArg, formatArg, drlftArg, endArg1};
    if (arg_nullcheck(normalArgTable) != 0)
    {
        printf("error allocating normalArgTable\n");
        return ERROR;
    }

#if defined(_WIN32) || defined(_WIN64)
    struct arg_str  *cmdArg  = arg_str1 (NULL, NULL, "compare", NULL);
#else
    struct arg_rex  *cmdArg  = arg_rex1 (NULL, NULL, "compare", NULL, REG_ICASE, NULL);
#endif
    struct arg_file *compArg = arg_filen(NULL, NULL, NULL, 1, 2, "result file(s)");
    struct arg_end  *endArg2 = arg_end  (10);

    void *compareArgTable[] = {cmdArg, compArg, endArg2};
    if (arg_nullcheck(compareArgTable) != 0)
    {
        printf("error allocating compareArgTable\n");
        return ERROR;
    }

    ArgResult result = NORMAL;

    nodeArg->filename[0] = inFn->node;
    bcArg->filename[0]   = inFn->bc;
    initArg->filename[0] = inFn->init;
    outArg->filename[0]  = inFn->result;
    compArg->filename[1] = inFn->final;

    uArg->dval[0]   = args->u;
    vArg->dval[0]   = args->v;
    rhoArg->dval[0] = args->rho;
    visArg->dval[0] = args->viscosity;

    iterArg->ival[0]  = args->iterations;
    autoeArg->ival[0] = args->autosaveEvery;
    autoaArg->ival[0] = args->autosaveAfter;
    drlftArg->ival[0] = args->boundaryId;

    int normalErrorCount = arg_parse(argc, argv, normalArgTable);
    int compareErrorCount = arg_parse(argc, argv, compareArgTable);

    if (normalErrorCount == 0)
    {
        if (helpArg->count > 0)
        {
            printf("\nUsage: %s ", argv[0]);
            arg_print_syntax(stdout, normalArgTable, "\n");
            printf("Usage: %s ", argv[0]);
            arg_print_syntax(stdout, compareArgTable, "\n");
            printf("\n");
            arg_print_glossary(stdout, normalArgTable, "%-31s %s\n");
            result = HELP;
        }
        else if (testArg->count > 0)
        {
            result = TEST;
        }
        else if (initArg->count == 0)
        {
            strcpy(inFn->node  , nodeArg->filename[0]);
            strcpy(inFn->bc    , bcArg->filename[0]);
            strcpy(inFn->init  , initArg->filename[0]);
            strcpy(inFn->result, outArg->filename[0]);

            int l = strlen(inFn->result);
            if (l && inFn->result[l-1] != '/')
            {
                strcat(inFn->result, "/");
            }

            args->u         = uArg->dval[0];
            args->v         = vArg->dval[0];
            args->rho       = rhoArg->dval[0];
            args->viscosity = visArg->dval[0];

            args->iterations    = iterArg->ival[0];
            args->autosaveEvery = autoeArg->ival[0];
            args->autosaveAfter = autoaArg->ival[0];
            args->boundaryId    = drlftArg->ival[0];

            args->inletProfile   = getInletProfile  (inletArg->sval[0],  args->inletProfile);
            args->outletProfile  = getOutletProfile (outlArg->sval[0],   args->outletProfile);
            args->collisionModel = getCollisionModel(collArg->sval[0],   args->collisionModel);
            args->boundaryType   = getBoundaryType  (curvedArg->count,   args->boundaryType);
            args->outputFormat   = getOutputFormat  (formatArg->sval[0], args->outputFormat);
        }
        else
        {
            strcpy(inFn->node  , nodeArg->filename[0]);
            strcpy(inFn->bc    , bcArg->filename[0]);
            strcpy(inFn->init  , initArg->filename[0]);
            strcpy(inFn->result, outArg->filename[0]);
            result = INIT;
        }
    }
    else if (compareErrorCount == 0)
    {
        if (cmdArg->count == 1 && compArg->count > 0)
        {
            strcpy(inFn->comp  , compArg->filename[0]);
            strcpy(inFn->final , compArg->filename[1]);
            result = COMPARE;
        }
    }
    else //none of the arguments worked
    {
        if (cmdArg->count > 0)
        {
          printf("2. ");
          arg_print_errors(stdout, endArg2, argv[0]);
        }
        else
        {
          printf("1. ");
          arg_print_errors(stdout, endArg1, argv[0]);
        }
        printf("\nUsage: %s ", argv[0]);
        arg_print_syntax(stdout, normalArgTable, "\n");
        printf("Usage: %s ", argv[0]);
        arg_print_syntax(stdout, compareArgTable, "\n");
        result = ERROR;
    }

    arg_freetable(normalArgTable,sizeof(normalArgTable)/sizeof(normalArgTable[0]));
    arg_freetable(compareArgTable,sizeof(compareArgTable)/sizeof(compareArgTable[0]));
    return result;
}
