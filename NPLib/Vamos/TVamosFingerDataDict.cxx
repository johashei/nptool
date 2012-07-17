//
// File generated by rootcint at Tue Jul 17 16:41:37 2012

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME TVamosFingerDataDict
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "TVamosFingerDataDict.h"

#include "TClass.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"

// START OF SHADOWS

namespace ROOT {
   namespace Shadow {
   } // of namespace Shadow
} // of namespace ROOT
// END OF SHADOWS

namespace ROOT {
   void TVamosFingerData_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_TVamosFingerData(void *p = 0);
   static void *newArray_TVamosFingerData(Long_t size, void *p);
   static void delete_TVamosFingerData(void *p);
   static void deleteArray_TVamosFingerData(void *p);
   static void destruct_TVamosFingerData(void *p);
   static void streamer_TVamosFingerData(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TVamosFingerData*)
   {
      ::TVamosFingerData *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TVamosFingerData >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TVamosFingerData", ::TVamosFingerData::Class_Version(), "./TVamosFingerData.h", 31,
                  typeid(::TVamosFingerData), DefineBehavior(ptr, ptr),
                  &::TVamosFingerData::Dictionary, isa_proxy, 0,
                  sizeof(::TVamosFingerData) );
      instance.SetNew(&new_TVamosFingerData);
      instance.SetNewArray(&newArray_TVamosFingerData);
      instance.SetDelete(&delete_TVamosFingerData);
      instance.SetDeleteArray(&deleteArray_TVamosFingerData);
      instance.SetDestructor(&destruct_TVamosFingerData);
      instance.SetStreamerFunc(&streamer_TVamosFingerData);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TVamosFingerData*)
   {
      return GenerateInitInstanceLocal((::TVamosFingerData*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TVamosFingerData*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
TClass *TVamosFingerData::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *TVamosFingerData::Class_Name()
{
   return "TVamosFingerData";
}

//______________________________________________________________________________
const char *TVamosFingerData::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TVamosFingerData*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TVamosFingerData::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TVamosFingerData*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void TVamosFingerData::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TVamosFingerData*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *TVamosFingerData::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TVamosFingerData*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
void TVamosFingerData::Streamer(TBuffer &R__b)
{
   // Stream an object of class TVamosFingerData.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> fVamos_Finger_Energy;
      R__b.CheckByteCount(R__s, R__c, TVamosFingerData::IsA());
   } else {
      R__c = R__b.WriteVersion(TVamosFingerData::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b << fVamos_Finger_Energy;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

//______________________________________________________________________________
void TVamosFingerData::ShowMembers(TMemberInspector &R__insp)
{
      // Inspect the data members of an object of class TVamosFingerData.
      TClass *R__cl = ::TVamosFingerData::IsA();
      if (R__cl || R__insp.IsA()) { }
      R__insp.Inspect(R__cl, R__insp.GetParent(), "fVamos_Finger_Energy", &fVamos_Finger_Energy);
      TObject::ShowMembers(R__insp);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TVamosFingerData(void *p) {
      return  p ? new(p) ::TVamosFingerData : new ::TVamosFingerData;
   }
   static void *newArray_TVamosFingerData(Long_t nElements, void *p) {
      return p ? new(p) ::TVamosFingerData[nElements] : new ::TVamosFingerData[nElements];
   }
   // Wrapper around operator delete
   static void delete_TVamosFingerData(void *p) {
      delete ((::TVamosFingerData*)p);
   }
   static void deleteArray_TVamosFingerData(void *p) {
      delete [] ((::TVamosFingerData*)p);
   }
   static void destruct_TVamosFingerData(void *p) {
      typedef ::TVamosFingerData current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TVamosFingerData(TBuffer &buf, void *obj) {
      ((::TVamosFingerData*)obj)->::TVamosFingerData::Streamer(buf);
   }
} // end of namespace ROOT for class ::TVamosFingerData

/********************************************************
* TVamosFingerDataDict.cxx
* CAUTION: DON'T CHANGE THIS FILE. THIS FILE IS AUTOMATICALLY GENERATED
*          FROM HEADER FILES LISTED IN G__setup_cpp_environmentXXX().
*          CHANGE THOSE HEADER FILES AND REGENERATE THIS FILE.
********************************************************/

#ifdef G__MEMTEST
#undef malloc
#undef free
#endif

#if defined(__GNUC__) && __GNUC__ >= 4 && ((__GNUC_MINOR__ == 2 && __GNUC_PATCHLEVEL__ >= 1) || (__GNUC_MINOR__ >= 3))
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

extern "C" void G__cpp_reset_tagtableTVamosFingerDataDict();

extern "C" void G__set_cpp_environmentTVamosFingerDataDict() {
  G__add_compiledheader("TObject.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("TVamosFingerData.h");
  G__cpp_reset_tagtableTVamosFingerDataDict();
}
#include <new>
extern "C" int G__cpp_dllrevTVamosFingerDataDict() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* TVamosFingerData */
static int G__TVamosFingerDataDict_162_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   TVamosFingerData* p = NULL;
   char* gvp = (char*) G__getgvp();
   int n = G__getaryconstruct();
   if (n) {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new TVamosFingerData[n];
     } else {
       p = new((void*) gvp) TVamosFingerData[n];
     }
   } else {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new TVamosFingerData;
     } else {
       p = new((void*) gvp) TVamosFingerData;
     }
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__TVamosFingerDataDictLN_TVamosFingerData));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TVamosFingerDataDict_162_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((TVamosFingerData*) G__getstructoffset())->Clear();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TVamosFingerDataDict_162_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 114, (long) ((TVamosFingerData*) G__getstructoffset())->GetFingerEnergy());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TVamosFingerDataDict_162_0_6(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((TVamosFingerData*) G__getstructoffset())->SetFingerEnergy((UShort_t) G__int(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TVamosFingerDataDict_162_0_7(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) TVamosFingerData::Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TVamosFingerDataDict_162_0_8(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) TVamosFingerData::Class_Name());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TVamosFingerDataDict_162_0_9(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 115, (long) TVamosFingerData::Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TVamosFingerDataDict_162_0_10(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      TVamosFingerData::Dictionary();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TVamosFingerDataDict_162_0_14(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((TVamosFingerData*) G__getstructoffset())->StreamerNVirtual(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TVamosFingerDataDict_162_0_15(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) TVamosFingerData::DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TVamosFingerDataDict_162_0_16(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) TVamosFingerData::ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TVamosFingerDataDict_162_0_17(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) TVamosFingerData::ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TVamosFingerDataDict_162_0_18(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) TVamosFingerData::DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__TVamosFingerDataDict_162_0_19(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)

{
   TVamosFingerData* p;
   void* tmp = (void*) G__int(libp->para[0]);
   p = new TVamosFingerData(*(TVamosFingerData*) tmp);
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__TVamosFingerDataDictLN_TVamosFingerData));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef TVamosFingerData G__TTVamosFingerData;
static int G__TVamosFingerDataDict_162_0_20(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   char* gvp = (char*) G__getgvp();
   long soff = G__getstructoffset();
   int n = G__getaryconstruct();
   //
   //has_a_delete: 1
   //has_own_delete1arg: 0
   //has_own_delete2arg: 0
   //
   if (!soff) {
     return(1);
   }
   if (n) {
     if (gvp == (char*)G__PVOID) {
       delete[] (TVamosFingerData*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((TVamosFingerData*) (soff+(sizeof(TVamosFingerData)*i)))->~G__TTVamosFingerData();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (TVamosFingerData*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((TVamosFingerData*) (soff))->~G__TTVamosFingerData();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic assignment operator
static int G__TVamosFingerDataDict_162_0_21(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   TVamosFingerData* dest = (TVamosFingerData*) G__getstructoffset();
   *dest = *(TVamosFingerData*) libp->para[0].ref;
   const TVamosFingerData& obj = *dest;
   result7->ref = (long) (&obj);
   result7->obj.i = (long) (&obj);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* TVamosFingerData */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncTVamosFingerDataDict {
 public:
  G__Sizep2memfuncTVamosFingerDataDict(): p(&G__Sizep2memfuncTVamosFingerDataDict::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncTVamosFingerDataDict::*p)();
};

size_t G__get_sizep2memfuncTVamosFingerDataDict()
{
  G__Sizep2memfuncTVamosFingerDataDict a;
  G__setsizep2memfunc((int)a.sizep2memfunc());
  return((size_t)a.sizep2memfunc());
}


/*********************************************************
* virtual base class offset calculation interface
*********************************************************/

   /* Setting up class inheritance */

/*********************************************************
* Inheritance information setup/
*********************************************************/
extern "C" void G__cpp_setup_inheritanceTVamosFingerDataDict() {

   /* Setting up class inheritance */
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__TVamosFingerDataDictLN_TVamosFingerData))) {
     TVamosFingerData *G__Lderived;
     G__Lderived=(TVamosFingerData*)0x1000;
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__TVamosFingerDataDictLN_TVamosFingerData),G__get_linked_tagnum(&G__TVamosFingerDataDictLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
   }
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableTVamosFingerDataDict() {

   /* Setting up typedef entry */
   G__search_typename2("UShort_t",114,-1,0,-1);
   G__setnewtype(-1,"Unsigned Short integer 2 bytes (unsigned short)",0);
   G__search_typename2("Version_t",115,-1,0,-1);
   G__setnewtype(-1,"Class version identifier (short)",0);
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__TVamosFingerDataDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__TVamosFingerDataDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__TVamosFingerDataDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__TVamosFingerDataDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__TVamosFingerDataDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__TVamosFingerDataDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__TVamosFingerDataDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__TVamosFingerDataDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__TVamosFingerDataDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__TVamosFingerDataDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* TVamosFingerData */
static void G__setup_memvarTVamosFingerData(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__TVamosFingerDataDictLN_TVamosFingerData));
   { TVamosFingerData *p; p=(TVamosFingerData*)0x1000; if (p) { }
   G__memvar_setup((void*)0,114,0,0,-1,G__defined_typename("UShort_t"),-1,4,"fVamos_Finger_Energy=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__TVamosFingerDataDictLN_TClass),-1,-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarTVamosFingerDataDict() {
}
/***********************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
***********************************************************/

/*********************************************************
* Member function information setup for each class
*********************************************************/
static void G__setup_memfuncTVamosFingerData(void) {
   /* TVamosFingerData */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__TVamosFingerDataDictLN_TVamosFingerData));
   G__memfunc_setup("TVamosFingerData",1583,G__TVamosFingerDataDict_162_0_1, 105, G__get_linked_tagnum(&G__TVamosFingerDataDictLN_TVamosFingerData), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Clear",487,G__TVamosFingerDataDict_162_0_2, 121, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Clear",487,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "C - 'Option_t' 10 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Dump",406,(G__InterfaceMethod) NULL,121, -1, -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("GetFingerEnergy",1509,G__TVamosFingerDataDict_162_0_5, 114, -1, G__defined_typename("UShort_t"), 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("SetFingerEnergy",1521,G__TVamosFingerDataDict_162_0_6, 121, -1, -1, 0, 1, 1, 1, 0, "r - 'UShort_t' 0 - E", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Class",502,G__TVamosFingerDataDict_162_0_7, 85, G__get_linked_tagnum(&G__TVamosFingerDataDictLN_TClass), -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (TClass* (*)())(&TVamosFingerData::Class) ), 0);
   G__memfunc_setup("Class_Name",982,G__TVamosFingerDataDict_162_0_8, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&TVamosFingerData::Class_Name) ), 0);
   G__memfunc_setup("Class_Version",1339,G__TVamosFingerDataDict_162_0_9, 115, -1, G__defined_typename("Version_t"), 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (Version_t (*)())(&TVamosFingerData::Class_Version) ), 0);
   G__memfunc_setup("Dictionary",1046,G__TVamosFingerDataDict_162_0_10, 121, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (void (*)())(&TVamosFingerData::Dictionary) ), 0);
   G__memfunc_setup("IsA",253,(G__InterfaceMethod) NULL,85, G__get_linked_tagnum(&G__TVamosFingerDataDictLN_TClass), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ShowMembers",1132,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TMemberInspector' - 1 - insp", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Streamer",835,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - b", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("StreamerNVirtual",1656,G__TVamosFingerDataDict_162_0_14, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - b", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DeclFileName",1145,G__TVamosFingerDataDict_162_0_15, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&TVamosFingerData::DeclFileName) ), 0);
   G__memfunc_setup("ImplFileLine",1178,G__TVamosFingerDataDict_162_0_16, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&TVamosFingerData::ImplFileLine) ), 0);
   G__memfunc_setup("ImplFileName",1171,G__TVamosFingerDataDict_162_0_17, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&TVamosFingerData::ImplFileName) ), 0);
   G__memfunc_setup("DeclFileLine",1152,G__TVamosFingerDataDict_162_0_18, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&TVamosFingerData::DeclFileLine) ), 0);
   // automatic copy constructor
   G__memfunc_setup("TVamosFingerData", 1583, G__TVamosFingerDataDict_162_0_19, (int) ('i'), G__get_linked_tagnum(&G__TVamosFingerDataDictLN_TVamosFingerData), -1, 0, 1, 1, 1, 0, "u 'TVamosFingerData' - 11 - -", (char*) NULL, (void*) NULL, 0);
   // automatic destructor
   G__memfunc_setup("~TVamosFingerData", 1709, G__TVamosFingerDataDict_162_0_20, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 1);
   // automatic assignment operator
   G__memfunc_setup("operator=", 937, G__TVamosFingerDataDict_162_0_21, (int) ('u'), G__get_linked_tagnum(&G__TVamosFingerDataDictLN_TVamosFingerData), -1, 1, 1, 1, 1, 0, "u 'TVamosFingerData' - 11 - -", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncTVamosFingerDataDict() {
}

/*********************************************************
* Global variable information setup for each class
*********************************************************/
static void G__cpp_setup_global0() {

   /* Setting up global variables */
   G__resetplocal();

}

static void G__cpp_setup_global1() {

   G__resetglobalenv();
}
extern "C" void G__cpp_setup_globalTVamosFingerDataDict() {
  G__cpp_setup_global0();
  G__cpp_setup_global1();
}

/*********************************************************
* Global function information setup for each class
*********************************************************/
static void G__cpp_setup_func0() {
   G__lastifuncposition();

}

static void G__cpp_setup_func1() {
}

static void G__cpp_setup_func2() {
}

static void G__cpp_setup_func3() {
}

static void G__cpp_setup_func4() {
}

static void G__cpp_setup_func5() {
}

static void G__cpp_setup_func6() {
}

static void G__cpp_setup_func7() {
}

static void G__cpp_setup_func8() {
}

static void G__cpp_setup_func9() {
}

static void G__cpp_setup_func10() {
}

static void G__cpp_setup_func11() {
}

static void G__cpp_setup_func12() {

   G__resetifuncposition();
}

extern "C" void G__cpp_setup_funcTVamosFingerDataDict() {
  G__cpp_setup_func0();
  G__cpp_setup_func1();
  G__cpp_setup_func2();
  G__cpp_setup_func3();
  G__cpp_setup_func4();
  G__cpp_setup_func5();
  G__cpp_setup_func6();
  G__cpp_setup_func7();
  G__cpp_setup_func8();
  G__cpp_setup_func9();
  G__cpp_setup_func10();
  G__cpp_setup_func11();
  G__cpp_setup_func12();
}

/*********************************************************
* Class,struct,union,enum tag information setup
*********************************************************/
/* Setup class/struct taginfo */
G__linked_taginfo G__TVamosFingerDataDictLN_TClass = { "TClass" , 99 , -1 };
G__linked_taginfo G__TVamosFingerDataDictLN_TBuffer = { "TBuffer" , 99 , -1 };
G__linked_taginfo G__TVamosFingerDataDictLN_TMemberInspector = { "TMemberInspector" , 99 , -1 };
G__linked_taginfo G__TVamosFingerDataDictLN_TObject = { "TObject" , 99 , -1 };
G__linked_taginfo G__TVamosFingerDataDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__TVamosFingerDataDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__TVamosFingerDataDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__TVamosFingerDataDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__TVamosFingerDataDictLN_TVamosFingerData = { "TVamosFingerData" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableTVamosFingerDataDict() {
  G__TVamosFingerDataDictLN_TClass.tagnum = -1 ;
  G__TVamosFingerDataDictLN_TBuffer.tagnum = -1 ;
  G__TVamosFingerDataDictLN_TMemberInspector.tagnum = -1 ;
  G__TVamosFingerDataDictLN_TObject.tagnum = -1 ;
  G__TVamosFingerDataDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__TVamosFingerDataDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__TVamosFingerDataDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__TVamosFingerDataDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__TVamosFingerDataDictLN_TVamosFingerData.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableTVamosFingerDataDict() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__TVamosFingerDataDictLN_TClass);
   G__get_linked_tagnum_fwd(&G__TVamosFingerDataDictLN_TBuffer);
   G__get_linked_tagnum_fwd(&G__TVamosFingerDataDictLN_TMemberInspector);
   G__get_linked_tagnum_fwd(&G__TVamosFingerDataDictLN_TObject);
   G__get_linked_tagnum_fwd(&G__TVamosFingerDataDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__TVamosFingerDataDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__TVamosFingerDataDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__TVamosFingerDataDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__TVamosFingerDataDictLN_TVamosFingerData),sizeof(TVamosFingerData),-1,29952,"VamosFingerData structure",G__setup_memvarTVamosFingerData,G__setup_memfuncTVamosFingerData);
}
extern "C" void G__cpp_setupTVamosFingerDataDict(void) {
  G__check_setup_version(30051515,"G__cpp_setupTVamosFingerDataDict()");
  G__set_cpp_environmentTVamosFingerDataDict();
  G__cpp_setup_tagtableTVamosFingerDataDict();

  G__cpp_setup_inheritanceTVamosFingerDataDict();

  G__cpp_setup_typetableTVamosFingerDataDict();

  G__cpp_setup_memvarTVamosFingerDataDict();

  G__cpp_setup_memfuncTVamosFingerDataDict();
  G__cpp_setup_globalTVamosFingerDataDict();
  G__cpp_setup_funcTVamosFingerDataDict();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncTVamosFingerDataDict();
  return;
}
class G__cpp_setup_initTVamosFingerDataDict {
  public:
    G__cpp_setup_initTVamosFingerDataDict() { G__add_setup_func("TVamosFingerDataDict",(G__incsetup)(&G__cpp_setupTVamosFingerDataDict)); G__call_setup_funcs(); }
   ~G__cpp_setup_initTVamosFingerDataDict() { G__remove_setup_func("TVamosFingerDataDict"); }
};
G__cpp_setup_initTVamosFingerDataDict G__cpp_setup_initializerTVamosFingerDataDict;
