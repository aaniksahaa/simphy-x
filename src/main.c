/* SPDX-License-Identifier: GPL-3.0-only */
#include "sx_cli.h"
#include "sx_run.h"
#include "sx_tpool.h"
#include "sx_port.h"

int main(int argc, char **argv)
{
    sx_port_init();
    sx_limits_from_env();
    if (NUM_BUFFER < 20 || MAX_LEAVES < 100 || MAX_NAME < 10 || MAX_IT < 10 || IO_BUFFER < 100) {
        fprintf(stderr, "Environment variables have been modified to unproper values, please, do not change them unless it is completely necessary\n");
        return EXIT_FAILURE;
    }
    sx_settings *s = sx_settings_new();
    long prc = sx_get_settings(s, argc, argv);
    if (prc == SX_TERMINATE_NOERROR) { fflush(stdout); return EXIT_SUCCESS; }
    sx_check(prc, ": getting settings");
    int nthreads = s->threads > 0 ? s->threads : sx_hardware_threads();
    printf("%s", "   _____ _           ____  __             __  __\n  / ___/(_)___ ___  / __ \\/ /_  __  __    \\ \\/ /\n  \\__ \\/ / __ `__ \\/ /_/ / __ \\/ / / /_____\\  / \n ___/ / / / / / / / ____/ / / / /_/ /_____/  \\ \n/____/_/_/ /_/ /_/_/   /_/ /_/\\__, /     /_/\\_\\\n                             /____/\n");
    printf("  simphy-x %s  |  mode: %s  |  %d thread%s\n\n---------------------------------------------------\n\n", SX_VERSION, s->mode == SX_MODE_FAST ? "fast (independent random streams per locus)" : "compat (byte-identical to SimPhy 1.0.2)", nthreads, nthreads == 1 ? "" : "s");
    printf("\nGetting settings from command line...");
    printf("Done \n");
    long rc = sx_run(s);
    sx_check(rc, NULL);
    sx_settings_free(s);
    return EXIT_SUCCESS;
}
