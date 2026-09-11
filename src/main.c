#include "sx_cli.h"
#include "sx_run.h"

int main(int argc, char **argv)
{
    sx_limits_from_env();
    if (NUM_BUFFER < 20 || MAX_LEAVES < 100 || MAX_NAME < 10 || MAX_IT < 10 || IO_BUFFER < 100) {
        fprintf(stderr, "Environment variables have been modified to unproper values, please, do not change them unless it is completely necessary\n");
        return EXIT_FAILURE;
    }
    printf("%s", "   _____ _           ____  __             __  __\n  / ___/(_)___ ___  / __ \\/ /_  __  __    \\ \\/ /\n  \\__ \\/ / __ `__ \\/ /_/ / __ \\/ / / /_____\\  / \n ___/ / / / / / / / ____/ / / / /_/ /_____/  \\ \n/____/_/_/ /_/ /_/_/   /_/ /_/\\__, /     /_/\\_\\\n                             /____/  simphy-x " SX_VERSION "\n\n---------------------------------------------------\n\n");
    printf("\nGetting settings from command line...");
    sx_settings *s = sx_settings_new();
    sx_check(sx_get_settings(s, argc, argv), ": getting settings");
    printf("Done \n");
    long rc = sx_run(s);
    sx_check(rc, NULL);
    sx_settings_free(s);
    return EXIT_SUCCESS;
}
