// directories.h
// create directories
#pragma once

#ifdef _WIN64
#include <direct.h>
#include <sys/stat.h>
#define STAT _stat
#define MKDIR _mkdir

#else
#include <sys/stat.h>
#define MKDIR mkdir
#define STAT stat
#endif
#include <iostream>
// dir name
#include <string>
// MPI_RANK(x)
#include "parallel/mpi_manager.h"

struct DirCreator {
  // use: DirCreator::Create_Dir("./dir"); or use absolute path
  // this will create DirName in the same directory as the executable file
  // DO NOT FORGET the "." at the beginning of the path if relative path is used
  // DO NOT create directory with sub-directory in one call, e.g. "./dir/subdir"
  static void Create_Dir(const std::string& DirName) {
    MPI_RANK(0)
    struct STAT info;
    if (STAT(DirName.c_str(), &info) != 0) {
      // win: _mkdir; linux: mkdir
#ifdef _WIN64
      int status = MKDIR(DirName.c_str());
#else
      int status = MKDIR(DirName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif
      if (status == -1) {
        std::cout << "Error Creating Directory: " << DirName << std::endl;
      } else {
        std::cout << "Created Directory: " << DirName << std::endl;
      }
    }
  }
  // variadic template is used to create multiple directories
  // use: DirCreator::Create_Dir("./dir", "./dir/subdir1", "./dir/subdir2"...);
  // this will create DirName in the same directory as the executable file
  // DO NOT forget the "." at the beginning of the path if relative path is used
  // this function is called recursively:
  // e.g., DirCreator::Create_Dir("./dir", "./dir1", "./dir2");
  // 1. DirName = "./dir", args = {"./dir1", "./dir2"}
  //    Create_Dir("./dir"); -> Create_Dir(std::string DirName)
  //    Create_Dir("./dir1", "./dir2"); -> Create_Dir(std::string DirName,
  //    Args... args)
  // 2. DirName = "./dir1", args = {"./dir2"}
  //    Create_Dir("./dir1"); -> Create_Dir(std::string DirName)
  //    Create_Dir("./dir2"); -> Create_Dir(std::string DirName)
  template <typename... Args>
  static void Create_Dir(const std::string& DirName, const Args&... args) {
    Create_Dir(DirName);
    Create_Dir(args...);
  }
};

#include <filesystem>
struct DirCreator_ {
  static void Create_Dir(const std::string& DirName) {
    // Validate input
    if (DirName.empty()) {
      std::cerr << "Error: Directory name cannot be empty." << std::endl;
      return;
    }
    // if (DirName.find("..") != std::string::npos) {
    //   std::cerr << "Error: Directory name cannot contain '..'." << std::endl;
    //   return;
    // }
    // if (DirName.find_first_of("/\\:*?\"<>|") != std::string::npos) {
    //   std::cerr << "Error: Directory name contains invalid characters." <<
    //   std::endl; return;
    // }

    // Create directory
    try {
      std::filesystem::create_directory(DirName);
      std::cout << "Created Directory: " << DirName << std::endl;
    } catch (const std::filesystem::filesystem_error& e) {
      std::cerr << "Error Creating Directory: " << e.what() << std::endl;
    }
  }

  template <typename... Args>
  static void Create_Dir(const std::string& DirName, const Args&... args) {
    Create_Dir(DirName);
    Create_Dir(args...);
  }
};