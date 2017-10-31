#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdarg.h>
#include <stddef.h>
#include <unistd.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <time.h>
#include <limits.h>
#include <assert.h>

#if 0
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
#endif

#include "model.h"
#include "util_sdl.h"
#include "util_misc.h"

static void help(void);

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
    // XXX display_handler();

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

#if 0
//#define DEFAULT_WIN_WIDTH    1920
//#define DEFAULT_WIN_HEIGHT   1000
#define DEFAULT_WIN_WIDTH    1900
#define DEFAULT_WIN_HEIGHT   1000

#define FUSOR_PANE_SIZE       804

// XXX
char about[] = "\
000\n\
hello\n\
hello\n\
hello\n\
hello\n\
hello\n\
hello\n\
hello\n\
hello\n\
hello\n\
hello\n\
hello\n\
hello\n\
hello\n\
hello\n\
hello\n\
hello\n\
hello\n\
hello\n\
hello\n\
hello\n\
hello\n\
hello\n\
hello\n\
hello\n\
hello\n\
hello\n\
hello\n\
hello\n\
hello\n\
hello\n\
hello\n\
hello\n\
hello\n\
hello\n\
hello\n\
hello\n\
hello\n\
hello\n\
hello\n\
hello\n\
hello\n\
123";
static void draw_fusor_pane(rect_t * border_pane, rect_t * pane);

int32_t test_counter;

static void display_handler(void)
{
    //INFO("starting\n");
    //model_start();
    //pause();
    //INFO("done\n");

    int32_t       win_width, win_height;
    rect_t        fusor_pane_full, fusor_pane;
    bool          redraw, done=false;
    sdl_event_t * event;

    // XXX maybe sdl_init shoud draw the screen and process any events,and
    //     return the actual size
    win_width  = DEFAULT_WIN_WIDTH;
    win_height = DEFAULT_WIN_HEIGHT;
    if (sdl_init(&win_width, &win_height, true, false, NULL) < 0) {
        FATAL("sdl_init %dx%d failed\n", win_width, win_height);
    }

    DEBUG("XXX W H  %d %d\n", win_width, win_height);

    sdl_init_pane(&fusor_pane_full, &fusor_pane,
                  100, 100, 
                  FUSOR_PANE_SIZE, FUSOR_PANE_SIZE);

    // loop until done
    while (!done) {
        // XXX should the panes be inited after this
        sdl_get_state(&win_width, &win_height, NULL);
        DEBUG("win_width = %d win_height = %d\n", win_width, win_height);

        // initialize for display update
        sdl_display_init();

        // draw panes
        draw_fusor_pane(&fusor_pane_full, &fusor_pane);

        // register for events   
        sdl_event_register(SDL_EVENT_KEY_SHIFT_ESC, SDL_EVENT_TYPE_KEY, NULL);   // done (shift-esc key)
        sdl_event_register('?', SDL_EVENT_TYPE_KEY, NULL);
        sdl_event_register('r', SDL_EVENT_TYPE_KEY, NULL);
        sdl_event_register('q', SDL_EVENT_TYPE_KEY, NULL);

        // present the display
        sdl_display_present();

        // process events until either redraw or done flag is set
        redraw = false;
        while (true) {
            event = sdl_poll_event();
            if (event->event != SDL_EVENT_NONE) {
                DEBUG("XXX GOT EVENT %d\n", event->event);
            }
            switch (event->event) {
            case SDL_EVENT_QUIT: case SDL_EVENT_KEY_SHIFT_ESC: case 'q':
                done = true;
                break;
            case '?':
                sdl_display_text(about);
                redraw = true;
                break;
            case 'r':
                redraw = true;
                break;
            case SDL_EVENT_SCREENSHOT_TAKEN:
                redraw = true;
                break;
            case SDL_EVENT_WIN_SIZE_CHANGE:
            case SDL_EVENT_WIN_RESTORED:
                redraw = true;
                break;
            case SDL_EVENT_USER_START:
                test_counter++;
                redraw = true;
                break;
            default:
                break;
            }

            if (done) {
                return;
            }
            if (redraw) {
                break;
            }
            usleep(1000);
        }
    }
}

static void draw_fusor_pane(rect_t * border_pane, rect_t * pane)
{
    // draw the border
    DEBUG("borderpane %d %d %d %d\n",
        border_pane->x,
        border_pane->y,
        border_pane->w,
        border_pane->h);
    sdl_render_pane_border(border_pane, GREEN);

    // pane
    DEBUG("pane %d %d %d %d\n",
        pane->x,
        pane->y,
        pane->w,
        pane->h);

    // draw circles to represent the chamber and the grid
    //sdl_render_text(pane, 0, 0, 0, "LIVE12345678", GREEN, BLACK);
    sdl_render_text_with_event(pane, 0, 0, 0, "LIVE1234567890_1234567890_1234567890_1234567890_1234567890_ABCDEFG", GREEN, BLACK, SDL_EVENT_USER_START);
    sdl_render_circle(pane, 400, 400, 425, 2, WHITE);

    char s[100];
    sprintf(s, "%d", test_counter);
    sdl_render_text(pane, 3, 0, 0, s, GREEN, BLACK);

    return;

    // XXX test
    texture_t t = sdl_create_filled_circle_texture(100, RED);
    sdl_render_texture(pane, t, 0, 0);
    sdl_render_texture(pane, t, 700, 700);

    rect_t dst = {300,250,200,300};
    sdl_render_scaled_texture(pane, t, &dst);
    
    rect_t dst1 = {300,700,200,300};
    sdl_render_scaled_texture(pane, t, &dst1);

    rect_t dst2 = {700,700,200,300};
    sdl_render_scaled_texture(pane, t, &dst2);

    rect_t dst3 = {700,-200,200,300};
    sdl_render_scaled_texture(pane, t, &dst3);

    rect_t dst4 = {-100,-200,200,300};
    sdl_render_scaled_texture(pane, t, &dst4);

    rect_t dst5 = {-100,700,200,300};
    sdl_render_scaled_texture(pane, t, &dst5);

    rect_t dst6 = {-1,100,1600,300};
    sdl_render_scaled_texture(pane, t, &dst6);

    rect_t dst7 = {0,500,1600,300};
    sdl_render_scaled_texture(pane, t, &dst7);

    sdl_destroy_texture(t);



    // draw the particles, use red for ion and green for atom;
    // use larger size when zoomed in
}

//draw_circle(rect_t * pane) 
//{
    //// on first call create table of sin values
//}

#endif
