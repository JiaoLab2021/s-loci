#ifndef shell_hpp
#define shell_hpp
#include <stdint.h>
#include <string>
#include <string.h>
#include <utility>
#include <vector>

using namespace std;

namespace shell{
    /*
     *Environment£º
     *Linux(Ubuntu), C++11, gcc 7.5.0, g++ 7.5.0
     *Description£º
     *Execute Linux shell command and get command return value or command execution result
    */

    #ifndef PARAMETER_FLOW
    #define PARAMETER_FLOW
    #define IN
    #define OUT
    #define INOUT
    #endif    //PARAMETER_FLOW

    #ifndef BASE_TYPE_DEF
    #define BASE_TYPE_DEF
    
    typedef int16_t          SHORT;
    typedef uint16_t         USHORT;
    typedef int32_t          INT;
    typedef uint32_t         UINT;
    typedef int64_t          DLONG;
    typedef uint64_t         DULONG;
    typedef void             VOID;
    typedef bool             BOOL;
    typedef char             CHAR;
    typedef unsigned char    UCHAR;
    typedef float            FLOAT;
    typedef double           DOUBLE;
    #endif    //BASE_TYPE_DEF

    using std::make_pair;
    using std::pair;
    using std::string;
    using std::vector;

    class Shell {
    public:

        /*
        *Function    : exeShellCmd
        *Description : Execute Linux shell commands and get command return values
        *Input       : IN const string& cmd = "", Linux shell command
        *            : OUT INT* cmdReturnValue = nullptr, command return value
        *Return      : pair<BOOL, string>, <whether the function is executed successfully, error message when the function fails>
        *Caution     :
        */
        static pair<BOOL, string> exeShellCmd(IN const string& cmd = "", OUT INT* cmdReturnValue = nullptr);

        /*
        *Function    : exeShellCmd
        *Description : Execute Linux shell commands and get command execution results
        *Input       : IN const string& cmd, Linux shell command
        *            : OUT vector<string>& results, command execution results
        *Return      : pair<BOOL, string>, <whether the function is executed successfully, error message when the function fails>
        *Caution     :
        */
        static pair<BOOL, string> exeShellCmd(IN const string& cmd, OUT vector<string>& results);

    }; //Shell

    pair<BOOL, string> exeShellCmd(IN const string &cmd, OUT INT *cmdReturnValue) {
        pid_t status;                 // pid_t is int
        status = system(cmd.c_str()); // Phase 1: Create a child process and other preparations. If it fails, return -1
        if (-1 == status) {
            return make_pair(false, "Error: stage 1: " + string(strerror(errno)));
        } else {
            // Phase 2: Call /bin/sh to start the script execution. If the script fails to start or the script does not end normally, the reason value is written to the lower 8~15 bits of status.
            // No matter what value is returned in the script, whether it is 0 or non-0, it is considered a normal execution. Even if the script does not exist or has no execution permission, it is considered a normal execution.
            // If the script is forcibly killed during execution, it is considered an abnormal end.
            if (WIFEXITED(status)) {
                if (nullptr != cmdReturnValue) {
                    *cmdReturnValue = WEXITSTATUS(status); // Get the script return value. Generally, if the script or command is executed correctly, the return value is 0. If it is executed incorrectly, other values ??are returned.
                }
                return make_pair(true, "");
            } else {
                return make_pair(false, "Error: stage 2.");
            }
        }
    } //exeShellCmd()

    pair<BOOL, string> exeShellCmd(IN const string &cmd, OUT vector<string> &results) {
        INT bufferSize = 10240;
        CHAR *buffer = new CHAR[bufferSize];
        FILE *pFile = NULL;
        if (NULL == (pFile = popen(cmd.c_str(), "r"))) {
            return make_pair(false, "Execute shell command error");
        }
        while (NULL != fgets(buffer, bufferSize, pFile)) {
            buffer[strlen(buffer) - 1] = '\0';  // fgets() will automatically add a newline character at the end. In Linux, the newline character is \n(LF). Here, remove the automatically added newline character.
            results.emplace_back(buffer);
        }
        delete[] buffer;
        pclose(pFile);
        return make_pair(true, "");
    } //exeShellCmd()
}

#endif