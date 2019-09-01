/*!
        @file    action_list.cpp

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "action_list.h"

const std::string ActionList::class_name = "ActionList";

//====================================================================
ActionList::ActionList(const int Nlevels)
  : m_vl(CommonParameters::Vlevel()),
  m_Nlevels(Nlevels), m_actions(Nlevels), m_integrator_types(Nlevels)
{
}


//====================================================================
bool ActionList::append(const int level, Action *action)
{
  if ((level < 0) || (level >= m_Nlevels)) {
    vout.crucial(m_vl, "Error at %s: (level=%d > 0) && (level=%d < Nlevels=%d) is required\n", class_name.c_str(), level, level, m_Nlevels);
    exit(EXIT_FAILURE);
  }

  m_actions[level].push_back(action);

  return true;
}


//====================================================================
bool ActionList::append(const int level, unique_ptr<Action>& action)
{
  return append(level, action.get());
}


//====================================================================
bool ActionList::append(const int level, const std::vector<Action *>& actions)
{
  if ((level < 0) || (level >= m_Nlevels)) {
    vout.crucial(m_vl, "Error at %s: (level=%d > 0) && (level=%d < Nlevels=%d) is required\n", class_name.c_str(), level, level, m_Nlevels);
    exit(EXIT_FAILURE);
  }

  for (size_t i = 0, n = actions.size(); i < n; ++i) {
    m_actions[level].push_back(actions[i]);
  }

  return true;
}


//====================================================================
bool ActionList::set_integrator_type(const int level, const std::string& type)
{
  if ((level < 0) || (level >= m_Nlevels)) {
    vout.crucial(m_vl, "Error at %s: (level=%d > 0) && (level=%d < Nlevels=%d) is required\n", class_name.c_str(), level, level, m_Nlevels);
    exit(EXIT_FAILURE);
  }

  m_integrator_types[level] = type;

  return true;
}


//====================================================================
bool ActionList::set_integrator_type(const std::string& type)
{
  for (size_t i = 0; i < m_Nlevels; ++i) {
    m_integrator_types[i] = type;
  }

  return true;
}


//====================================================================
ActionSet ActionList::get_actions() const
{
  ActionSet v;

  for (size_t i = 0; i < m_Nlevels; ++i) {
    v.insert(v.end(), m_actions[i].begin(), m_actions[i].end());
  }

  return v;
}


//====================================================================
ActionSet ActionList::get_actions(const int level) const
{
  if ((level < 0) || (level >= m_Nlevels)) return ActionSet();

  return m_actions[level];
}


//====================================================================
std::string ActionList::get_integrator_type(const int level) const
{
  if ((level < 0) || (level >= m_Nlevels)) return std::string();

  return m_integrator_types[level];
}


//====================================================================
int ActionList::get_levels() const
{
  return m_Nlevels;
}


//====================================================================
//============================================================END=====
