// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dIhomedInstardIths_HH4bdIlimitSettingdIstep1dIdrawAlpha2_C_ACLiC_dict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "/home/nstar/ths_HH4b/limitSetting/step1/./drawAlpha2.C"

// Header files passed via #pragma extra_include

namespace {
  void TriggerDictionaryInitialization_drawAlpha2_C_ACLiC_dict_Impl() {
    static const char* headers[] = {
"./drawAlpha2.C",
0
    };
    static const char* includePaths[] = {
"/home/nstar/root-6.08.06/include",
"/home/nstar/root-6.08.06/etc",
"/home/nstar/root-6.08.06/etc/cling",
"/home/nstar/root-6.08.06/include",
"/usr/include/c++/6",
"/usr/include/x86_64-linux-gnu/c++/6",
"/usr/include/c++/6/backward",
"/home/nstar/root-6.08.06/include",
"/home/nstar/ths_HH4b/limitSetting/step1/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "drawAlpha2_C_ACLiC_dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "drawAlpha2_C_ACLiC_dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif
#ifndef __ACLIC__
  #define __ACLIC__ 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "./drawAlpha2.C"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"drawAlpha2", payloadCode, "@",
"drawAlphaBase", payloadCode, "@",
"setTH", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("drawAlpha2_C_ACLiC_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_drawAlpha2_C_ACLiC_dict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_drawAlpha2_C_ACLiC_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_drawAlpha2_C_ACLiC_dict() {
  TriggerDictionaryInitialization_drawAlpha2_C_ACLiC_dict_Impl();
}
