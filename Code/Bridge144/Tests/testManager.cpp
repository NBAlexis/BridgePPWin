/*!
        @file    $Id: testManager.cpp #$

        @brief

        @author  Shinji Motoki  (smotoki)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#include "Tests/testManager.h"

const std::string TestManager::class_name = "TestManager";


// prototype definition
std::vector<std::string> string_tokenize(const std::string& src, const char delim = '.');

//====================================================================
int DoNothing()
{
  return int();
}


//====================================================================
void TestManager::banner()
{
  vout.general(m_vl, "\n");
  vout.general(m_vl, "-------------------------------------------------\n");
  vout.general(m_vl, "------------ <Bridge++ %s> Test Menu -----------\n", BRIDGE_VERSION);
  vout.general(m_vl, "-------------------------------------------------\n");
  vout.general(m_vl, "Please select test category name \n");
}


//====================================================================
TestManager& TestManager::Instance()
{
  if (!s_instance) {
    //
    // acquire lock here
    //
    if (!s_instance) {
      CreateInstance();
    }
    //
    // release lock here
    //
  }

  return *s_instance;
}


//====================================================================
void TestManager::CreateInstance()
{
  static TestManager instance_;

  s_instance = &instance_;
}


//====================================================================
TestManager::TestManager()
  : m_vl(CommonParameters::Vlevel()),
    m_root_node("<top_level>"),
    m_precision(Test::default_precision)
{
}


//====================================================================
TestManager::~TestManager()
{
  s_instance = 0;
}


//====================================================================
void TestManager::interactive()
{
  if (Communicator::is_primary()) {
    interactive_master();
  } else {
    interactive_slave();
  }
}


//====================================================================
void TestManager::interactive_master()
{
  menu(&m_root_node, true);

  std::string jobname = "<terminate>";
  Communicator::broadcast(1, jobname, 0);

  vout.detailed(m_vl, "TestManager::interactive_master: rank=%d, exit.\n", Communicator::self());
}


//====================================================================
void TestManager::interactive_slave()
{
  for ( ; ; ) {
    std::string jobname;
    Communicator::broadcast(1, jobname, 0);

    vout.detailed(m_vl, "TestManager::interactive_slave: rank=%d, job=%s\n", Communicator::self(), jobname.c_str());

    if (jobname == "<terminate>") break;

    batch(jobname);
  }

  vout.detailed(m_vl, "TestManager::interactive_slave: rank=%d, exit.\n", Communicator::self());
}


//====================================================================
void TestManager::batch(const std::string& arg)
{
  Node *p = find_node(&m_root_node, string_tokenize(arg));

  if (p && (p->m_function != DoNothing)) {
    int result = p->m_function();
  } else {
    vout.general(m_vl, "%s: No such test is found: %s\n\n", class_name.c_str(), arg.c_str());
  }
}


//====================================================================
void TestManager::batch_recursive(const std::string& arg)
{
  vout.detailed(m_vl, "TestManager::batch_recursive: arg=\"%s\"\n", arg.c_str());

  if (Communicator::is_primary()) {
    Node *p = (arg == "") ? &m_root_node : find_node(&m_root_node, string_tokenize(arg));

    if (p) {
      m_stat.reset();

      run_traversal(p);

      stat_report();
    } else {
      vout.general(m_vl, "%s: No such test is found: %s\n\n", class_name.c_str(), arg.c_str());
    }

    // send terminate signal to slave nodes
    std::string jobname = "<terminate>";
    Communicator::broadcast(1, jobname, 0);

    vout.detailed(m_vl, "TestManager::batch_recursive: rank=%d, exit.\n", Communicator::self());
  } else {
    interactive_slave();
  }
}


//====================================================================
void TestManager::batch_recursive(const int argc, char **argv)
{
  vout.detailed(m_vl, "TestManager::batch_recursive: argc=%d\n", argc);

  if (Communicator::is_primary()) {
    m_stat.reset();

    for (int i = 0; i < argc; ++i) {
      vout.detailed(m_vl, "TestManager::batch_recursive: argv[%d]=%s\n", i, argv[i]);

      std::string arg = std::string(argv[i]);

      Node *p = (arg == "") ? &m_root_node : find_node(&m_root_node, string_tokenize(arg));

      if (p) {
        run_traversal(p);
      } else {
        vout.general(m_vl, "%s: No such test is found: %s\n\n", class_name.c_str(), arg.c_str());
      }
    }

    stat_report();

    // send terminate signal to slave nodes
    std::string jobname = "<terminate>";
    Communicator::broadcast(1, jobname, 0);

    vout.detailed(m_vl, "TestManager::batch_recursive: rank=%d, exit.\n", Communicator::self());
  } else {
    interactive_slave();
  }
}


//====================================================================
bool TestManager::registerTest(const std::string& key, const Test_function func)
{
  append_key(&m_root_node, string_tokenize(key))->m_function = func;
  return true;
}


//====================================================================
void TestManager::run(const Node *p)
{
  if (!p) return;

  if (p->m_function == DoNothing) return;

  std::string testname = find_fullpath(p);
  vout.general(m_vl, "run test \"%s\"\n", testname.c_str());

  if (Communicator::is_primary()) {
    Communicator::broadcast(1, testname, 0);
  }

  int result = p->m_function();

  // increment stat counter
  if (result == -1) {
    m_stat.skip(testname);
  } else if (result == 0) {
    m_stat.success(testname);
  } else {
    m_stat.failure(testname);
  }

  // check result
  //check_result(testname, result);
}


//====================================================================
void TestManager::run_traversal(const Node *p)
{
  if (!p) return;

  if (p->m_next.size() == 0) {
    run(p);
  } else {
    for (size_t i = 0; i < p->m_next.size(); ++i) {
      run_traversal(p->m_next[i]);
    }
  }
}


//====================================================================
bool TestManager::menu(const Node *p, bool is_top)
{
  if (!p) return false;

  if (is_top && (p->m_next.size() > 0)) banner();


  if (p->m_next.size() == 0) {  // leaf node
    vout.general(m_vl, "run %s\n", p->m_name.c_str());
    run(p);
    return true;
  }

  bool do_continue = true;

  while (do_continue) // loop forever

  {                   // show item list
    for (unsigned int i = 0; i < p->m_next.size(); ++i) {
      vout.general(m_vl, "%u : %s\n", i + 1, p->m_next[i]->m_name.c_str());
    }

    vout.general(m_vl, "a : Test All \n");
    vout.general(m_vl, "p : Setup test check precision (current precision: %d)\n", m_precision);
    if (!is_top) {
      vout.general(m_vl, "u : Go back \n");
    }
    vout.general(m_vl, "q : Quit \n");

    char buf[1024];
    buf[0] = '\0';

    int choice = 0;

    bool do_alltest      = false;
    bool do_setprecision = false;

    while (true)    // loop until get valid answer.
    {
      vout.general(m_vl, "choice> ");
      // scanf("%1023s", buf);
      std::cin >> buf;

      if (buf[0] == 0) return false;    // ctrl-D to escape

      if (buf[0] == 'a') {
        do_alltest = true;
        break;
      }

      if (buf[0] == 'p') {
        do_setprecision = true;
        break;
      }

      if (!is_top) {
        if (buf[0] == 'u') return true;  // go up.
      }
      if (buf[0] == 'q') return false;   // quit.


      choice = atoi(buf);
      buf[0] = '\0';

      if ((choice >= 1) && (choice <= int(p->m_next.size()))) break;
    }

    if (do_alltest) {
      m_stat.reset();

      run_traversal(p);

      stat_report();
    } else if (do_setprecision) {
      set_precision();
    } else {
      do_continue = menu(p->m_next[choice - 1]);
    }

    if (!do_continue) break;
  }

  return do_continue;
}


//====================================================================
void TestManager::set_precision()
{
  vout.general(m_vl, "Please input test check precision \n");
  vout.general(m_vl, "input number from 1 to 14 (default: 12)\n");

  char buf[1024];
  buf[0] = '\0';

  vout.general(m_vl, "precision = ");
  // scanf("%1023s", buf);
  std::cin >> buf;

  int prec = atoi(buf);
  if ((prec >= 1) && (prec <= 14)) {
    m_precision = prec;

    vout.general(m_vl, "precision set to %d\n", m_precision);
  } else {
    vout.general(m_vl, "invalid value.\n");
  }
}


//====================================================================
TestManager::Node *TestManager::find_node(TestManager::Node *p, const std::vector<std::string>& v)
{
  if (!p) return 0;

  if (v.size() == 0) return p;

  for (std::vector<std::string>::const_iterator r = v.begin(); r != v.end(); ++r) {
    bool   is_found = false;
    size_t i        = 0;
    for (i = 0; i < p->m_next.size(); ++i) {
      if (p->m_next[i]->m_name == (*r)) {
        is_found = true;
        break;
      }
    }

    if (!is_found) return 0;

    p = p->m_next[i];
  }

  return p;
}


//====================================================================
TestManager::Node *TestManager::append_key(TestManager::Node *p, const std::vector<std::string>& v)
{
  if (!p) return 0;

  if (v.size() == 0) return p;

  for (std::vector<std::string>::const_iterator r = v.begin(); r != v.end(); ++r) {
    bool   is_found = false;
    size_t i        = 0;
    for (i = 0; i < p->m_next.size(); ++i) {
      if (p->m_next[i]->m_name == (*r)) {
        is_found = true;
        break;
      }
    }

    if (is_found) {
      p = p->m_next[i];
    } else {
      Node *q = new Node(*r, p);

      //append at the end of list
      // p->m_next.push_back(q);

      //append in alphabetical order
      std::vector<Node *>::iterator iter = p->m_next.begin();
      for ( ; iter != p->m_next.end(); ++iter) {
        if (*r < (*iter)->m_name) {
          break;
        }
      }
      p->m_next.insert(iter, q);

      p = q;
    }
  }

  return p;
}


//====================================================================
TestManager::Node *TestManager::append_key(TestManager::Node *p, int argc, char **argv)
{
  if (!p) return 0;

  if (argc == 0) return p;

  vout.general(m_vl, "%s: p=%p, argc=%d, argv={ ", __func__, p, argc);
  for (int i = 0; i < argc; ++i) {
    vout.general(m_vl, "%s, ", argv[i]);
  }
  vout.general(m_vl, "}\n");

  for (size_t i = 0; i < p->m_next.size(); ++i) {
    if (p->m_next[i]->m_name == argv[0]) {
      return append_key(p->m_next[i], argc - 1, argv + 1);
    }
  }

  Node *q = new Node(argv[0], p);
  p->m_next.push_back(q);

  return append_key(q, argc - 1, argv + 1);
}


//====================================================================
std::string TestManager::find_fullpath(const TestManager::Node *p, const std::string& path)
{
  if (!p) return path;

  if (!p->m_prev) return path;  // omit "<top_level>"

  if (path.length() > 0) {
    return find_fullpath(p->m_prev, p->m_name + '.' + path);
  } else {
    return find_fullpath(p->m_prev, p->m_name);
  }
}


//====================================================================
void TestManager::Stat::reset()
{
  m_num_tests   = 0;
  m_num_success = 0;
  m_num_failure = 0;
  m_num_skip    = 0;
  m_list_failure.resize(0);
  m_list_skip.resize(0);
}


//====================================================================
void TestManager::Stat::success(const std::string& test_name)
{
  ++m_num_tests;
  ++m_num_success;
}


//====================================================================
void TestManager::Stat::failure(const std::string& test_name)
{
  ++m_num_tests;
  ++m_num_failure;
  m_list_failure.push_back(test_name);
}


//====================================================================
void TestManager::Stat::skip(const std::string& test_name)
{
  ++m_num_tests;
  ++m_num_skip;
  m_list_skip.push_back(test_name);
}


//====================================================================
void TestManager::stat_report() const
{
  vout.general(m_vl, "Test report:\n");
  vout.general(m_vl, "  Total number of tests = %3d\n", m_stat.m_num_tests);
  vout.general(m_vl, "  Number of successes   = %3d\n", m_stat.m_num_success);
  vout.general(m_vl, "  Number of failures    = %3d\n", m_stat.m_num_failure);
  vout.general(m_vl, "  Number of skipped     = %3d\n", m_stat.m_num_skip);

  if (m_stat.m_num_failure > 0) {
    vout.general(m_vl, "  Failed tests:\n");
    for (std::vector<std::string>::const_iterator p = m_stat.m_list_failure.begin(); p != m_stat.m_list_failure.end(); ++p) {
      vout.general(m_vl, "    %s\n", p->c_str());
    }
  }

//   if (m_stat.m_num_skip > 0) {
//     vout.general(m_vl, "  Skipped tests:\n");
//     for (std::vector<std::string>::const_iterator p = m_stat.m_list_skip.begin(); p != m_stat.m_list_skip.end(); ++p) {
//       vout.general(m_vl, "    %s\n", p->c_str());
//     }
//   }

  if (m_stat.m_num_tests == 0) {
    vout.general(m_vl, "  No test is performed.\n");
  } else {
    if (m_stat.m_num_failure == 0) {
      vout.general(m_vl, "  All tests are performed successfully.\n");
    }
  }

  vout.general(m_vl, "\n");
}


//====================================================================
void TestManager::list_traverse(const Node *p, const string& prefix)
{
  if (!p) return;

  if (p->m_next.size() == 0) {
    if (prefix == "") {
      vout.general(m_vl, "%s\n", p->m_name.c_str());
    } else {
      vout.general(m_vl, "%s.%s\n", prefix.c_str(), p->m_name.c_str());
    }
    return;
  }

  string prefix_next = "";

  if (p == &m_root_node) {
  } else if (prefix == "") {
    prefix_next = p->m_name;
  } else {
    prefix_next = prefix + "." + p->m_name;
  }

  for (size_t i = 0; i < p->m_next.size(); ++i) {
    list_traverse(p->m_next[i], prefix_next);
  }
}


//====================================================================
void TestManager::list()
{
  list_traverse(&m_root_node, "");
}


//====================================================================
void TestManager::traverse(const Node *p, const std::string& indent)
{
  if (!p) return;

  if (p->m_next.size() == 0) {
    vout.general(m_vl, "%sleaf \"%s\"\n", indent.c_str(), p->m_name.c_str());
    return;
  }

  vout.general(m_vl, "%snode \"%s\"\n", indent.c_str(), p->m_name.c_str());

  for (size_t i = 0; i < p->m_next.size(); ++i) {
    traverse(p->m_next[i], indent + "  ");
  }
}


//====================================================================
std::vector<std::string> string_tokenize(const std::string& src, const char delim)
{
  std::vector<std::string> retv;

  size_t npos = src.length();
  size_t p    = 0;

  while (true)
  {
    size_t q = src.find(delim, p);

    if (q >= npos) {                 // not found
      retv.push_back(src.substr(p)); // append rest
      break;
    }

    retv.push_back(src.substr(p, q - p));
    p = q + 1;  // skip delimiter.
  }

  return retv;
}


//====================================================================
TestManager *TestManager::s_instance = 0;

//====================================================================
//==============================================================END===
