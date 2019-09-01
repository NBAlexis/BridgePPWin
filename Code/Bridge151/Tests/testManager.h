/*!
        @file    $Id: testManager.h #$

        @brief

        @author  Shinji Motoki  (smotoki)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#ifndef TESTMANAGER_INCLUDED
#define TESTMANAGER_INCLUDED

#include <map>

#include "Communicator/communicator.h"
#include "Parameters/commonParameters.h"
#include "Tests/test.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

// typedef Functor<int> Test_function;
typedef int (*Test_function)(void);

// null function
int DoNothing();


//! TestManager class for managing and performing tests.

/**
   TestManager class provides framework for managing and performing
   various tests. Tests are registered to TestManager, and performed
   interactively through hierarchical menu, or in a batch mode by
   specifying test names.

   Each test is in a form of a function which takes no argument,
   and returns status whether the test is successful or not.
   It is expected that each test is independent with each other,
   leaving no side effects.

   TestManager collects and reports statistics of success/failure of tests
   when 'all' of the subtree of tests are selected in interactive mode.

   N.B. TestManager is a singleton class.
 */

class TestManager
{
 public:
  static const std::string class_name;

 public:
  void interactive();
  void batch(const std::string& arg);
  void batch_recursive(const std::string& arg = "");
  void batch_recursive(const int argc, char **argv);

  bool registerTest(const std::string& key, const Test_function func);

  void list();

  static TestManager& Instance();

  static bool RegisterTest(const std::string& key, const Test_function func)
  { return Instance().registerTest(key, func); }

  // for test result verification
  int getCheckPrecision() const
  { return m_precision; }

 private:
  // singleton
  TestManager();

  TestManager(const TestManager&);            // not implemented
  TestManager& operator=(const TestManager&); // not implemented

  ~TestManager();

  static void CreateInstance();

  static TestManager *s_instance;

  // bi-directional tree representation of key
  struct Node
  {
    Node                *m_prev;
    std::vector<Node *> m_next;
    std::string         m_name;     // key
    Test_function       m_function; // payload

    Node(const std::string& name, Node *const prev = 0)
      : m_prev(prev), m_next(), m_name(name), m_function(DoNothing) {}
    ~Node()
    {
      for (size_t i = 0; i < m_next.size(); ++i) { delete m_next[i]; }
    }
  };

  Node *find_node(Node *p, const std::vector<std::string>& v);
  Node *append_key(Node *p, const std::vector<std::string>& v);
  Node *append_key(Node *p, int argc, char **argv);

  // statistics
  struct Stat
  {
    int                      m_num_tests;
    int                      m_num_success;
    int                      m_num_failure;
    int                      m_num_skip;
    std::vector<std::string> m_list_failure;
    std::vector<std::string> m_list_skip;

    void                     reset();
    void                     success(const std::string& test_name);
    void                     failure(const std::string& test_name);
    void                     skip(const std::string& test_name);
  };

  void stat_report() const;

  std::string find_fullpath(const Node *p, const std::string& path = "");

  void interactive_master();
  void interactive_slave();

  void run(const Node *p);
  void run_traversal(const Node *p);

  bool menu(const Node *p, const bool is_top = false);
  void banner();

  void set_precision();

  void list_traverse(const Node *p, const std::string& prefix);

  // for debug
  void traverse(const Node *p, const std::string& indent = "");


  Bridge::VerboseLevel m_vl;

  Node m_root_node;

  // for test result verification
  int m_precision;

  // for statistics report
  Stat m_stat;
};
#endif // TEST_MANAGER_H
