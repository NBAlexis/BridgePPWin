/*!
@file    $Id:: factory.h #$

@brief

@author  Tatsumi Aoyama  (aoym)
$LastChangedBy: aoym $

@date    $LastChangedDate:: 2017-02-24 00:49:46 #$

@version $LastChangedRevision: 1561 $
*/

#ifndef FACTORY_INCLUDED
#define FACTORY_INCLUDED

//! Factory template class

/**
FactoryTemplate class provides framework of factories for classes
which relate identifiers of subclasses (class name string) and
instance-creation callbacks (functions).

This template takes two template parameters: AbstractProduct represents
typename of abstract base class of the product class family, and
ProductCreator specifies the type of function (including types of
arguments).
*/

#include <string>
#include <map>

typedef std::string   IdentifierType;

template<class AbstractProduct, typename ProductCreator>
class FactoryTemplate
{
private:
    typedef std::map<IdentifierType, ProductCreator>   Map;


public:
    static ProductCreator Find(const IdentifierType& subtype)
    {
        return Instance().find(subtype);
    }

    static bool Register(const IdentifierType& subtype, ProductCreator creator)
    {
        return Instance().append(subtype, creator);
    }

    ProductCreator find(const IdentifierType& subtype)
    {
        typename Map::const_iterator i = m_map.find(subtype);
        if (i != m_map.end()) {
            return i->second;
        }
        else {
            fprintf(stderr, "Factory::find: unknown type \"%s\"\n", subtype.c_str());

            exit(EXIT_FAILURE);

            return 0;
        }
    }

    bool append(const IdentifierType& subtype, ProductCreator creator)
    {
        if ((m_map.find(subtype) == m_map.end()) &&
            m_map.insert(typename Map::value_type(subtype, creator)).second) {
            return true;
        }
        else {
            fprintf(stderr, "Factory::append: duplicate key \"%s\"\n", subtype.c_str());
            return false;
        }
    }

    static FactoryTemplate& Instance()
    {
        if (!s_instance) {
            // lock
            if (!s_instance) {
                create_instance();
            }
            // unlock
        }
        return *s_instance;
    }

private:
    Map m_map;

    // singleton
    FactoryTemplate() {}
    FactoryTemplate(const FactoryTemplate&) {}
    FactoryTemplate& operator=(const FactoryTemplate&);

    virtual ~FactoryTemplate() {}

    static inline void create_instance()
    {
        static FactoryTemplate instance;

        s_instance = &instance;
    }

    static FactoryTemplate *s_instance;
};

template<class AbstractProduct, typename ProductCreator>
FactoryTemplate<AbstractProduct, ProductCreator> *FactoryTemplate<AbstractProduct, ProductCreator>::s_instance = 0;
#endif
