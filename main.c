#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdarg.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <math.h> 
#include <assert.h>

#include "model.h"
#include "util_sdl.h"
#include "util_misc.h"

static void help(void);
static void display_handler(void);

// -----------------  MAIN  -------------------------------------------------------------

int main(int argc, char **argv)
{
    static char default_params_str[] = "150,38,15,30,5";
    char * params_str = NULL;
    char * filename_str = NULL;

    // get and process options
    // -p <params_str>   : example "150,38,15,30,5" 
    //                       150 = chamber diameter in mm
    //                       38  = grid diameter in mm
    //                       15  = pressure in mtorr
    //                       30  = voltage in kv
    //                       5   = current in ma
    // -f <filename_str> : saved model filename
    // -h                : help
    while (true) {
        static bool opt_p_supplied = false;
        static bool opt_f_supplied = false;
        char opt_char = getopt(argc, argv, "p:f:h");
        if (opt_char == -1) {
            break;
        }
        switch (opt_char) {
        case 'p':
            if (opt_p_supplied || opt_f_supplied) {
                FATAL("-p and -f can not be combined\n");
            }
            params_str = optarg;
            opt_p_supplied = true;
            break;
        case 'f':
            if (opt_p_supplied || opt_f_supplied) {
                FATAL("-p and -f can not be combined\n");
            }
            filename_str = optarg;
            opt_f_supplied = true;
            break;
        case 'h':
            help();
            exit(0);
        default:
            exit(1);
            break;
        }
    }

    // if neither params_str or filename supplied then
    // use default params_str
    if (!params_str && !filename_str) {
        params_str = default_params_str;
    }
 
    // initialize model
    if (params_str) {
        model_init_from_params(params_str);
    } else {
        model_init_from_file(filename_str);
    }

    // call display handler
    display_handler();

    // terminate the model
    model_terminate();

    // terminate 
    return 0;
}

static void help(void)
{
    // XXX tbd
    printf("TBD\n");
}

// -----------------  DISPLAY_HANDLER  --------------------------------------------------

static void display_handler(void)
{
    INFO("starting\n");
    model_start();
    pause();
    INFO("done\n");

}
