#pragma once

#if defined(SIMULATIONCPP_EXPORTS)
#  define DLLEXPORT __declspec(dllexport)
#else
#  define DLLEXPORT __declspec(dllimport)
#endif
