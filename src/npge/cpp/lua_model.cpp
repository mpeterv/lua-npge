/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cstdlib>
#include <memory>

#define LUA_LIB
#include <lua.hpp>

#include "model.hpp"
#include "throw_assert.hpp"

using namespace npge;

#define LUA_TRY try {

#define LUA_CATCH } catch (std::exception& e) { \
    return luaL_error(L, e.what()); \
} catch (...) { \
    return luaL_error(L, "Unknown error occurred"); \
}

// first upvalue: metatable for Sequence instance
LUALIB_API int lua_Sequence(lua_State *L) {
    size_t name_size, text_size;
    const char* name = luaL_checklstring(L, 1, &name_size);
    const char* text = luaL_checklstring(L, 2, &text_size);
    const char* description = "";
    size_t description_size = 0;
    if (lua_gettop(L) >= 3 && lua_type(L, 3) == LUA_TSTRING) {
        description = luaL_checklstring(L, 3,
                &description_size);
    }
    std::allocator<Sequence> a;
    std::auto_ptr<Sequence> seq(a.allocate(1));
LUA_TRY
    new (seq.get()) Sequence(name, description,
            text, text_size);
LUA_CATCH
    void* s = lua_newuserdata(L, sizeof(SequencePtr));
    new (s) SequencePtr(seq.release());
    // get metatable of Sequence
    lua_pushvalue(L, lua_upvalueindex(1));
    lua_setmetatable(L, -2);
    return 1;
}

static SequencePtr& lua_toseq(lua_State* L, int index) {
    void* v = luaL_checkudata(L, index, "npge_Sequence");
    SequencePtr* s = reinterpret_cast<SequencePtr*>(v);
    return *s;
}

LUALIB_API int lua_Sequence_gc(lua_State *L) {
    SequencePtr& seq = lua_toseq(L, 1);
    seq.reset();
    return 0;
}

LUALIB_API int lua_Sequence_type(lua_State *L) {
    lua_pushstring(L, "Sequence");
    return 1;
}

LUALIB_API int lua_Sequence_name(lua_State *L) {
    const SequencePtr& seq = lua_toseq(L, 1);
    const std::string& name = seq->name();
    lua_pushlstring(L, name.c_str(), name.size());
    return 1;
}

LUALIB_API int lua_Sequence_description(lua_State *L) {
    const SequencePtr& seq = lua_toseq(L, 1);
    const std::string& description = seq->description();
    lua_pushlstring(L, description.c_str(),
                    description.size());
    return 1;
}

LUALIB_API int lua_Sequence_genome(lua_State *L) {
    const SequencePtr& seq = lua_toseq(L, 1);
    std::string genome = seq->genome();
    if (!genome.empty()) {
        lua_pushlstring(L, genome.c_str(), genome.size());
    } else {
        lua_pushnil(L);
    }
    return 1;
}

LUALIB_API int lua_Sequence_chromosome(lua_State *L) {
    const SequencePtr& seq = lua_toseq(L, 1);
    std::string chr = seq->chromosome();
    if (!chr.empty()) {
        lua_pushlstring(L, chr.c_str(), chr.size());
    } else {
        lua_pushnil(L);
    }
    return 1;
}

LUALIB_API int lua_Sequence_circular(lua_State *L) {
    const SequencePtr& seq = lua_toseq(L, 1);
    int circular = seq->circular();
    lua_pushboolean(L, circular);
    return 1;
}

LUALIB_API int lua_Sequence_text(lua_State *L) {
    const SequencePtr& seq = lua_toseq(L, 1);
    const std::string& text = seq->text();
    lua_pushlstring(L, text.c_str(), text.size());
    return 1;
}

LUALIB_API int lua_Sequence_length(lua_State *L) {
    const SequencePtr& seq = lua_toseq(L, 1);
    lua_pushinteger(L, seq->length());
    return 1;
}

LUALIB_API int lua_Sequence_sub(lua_State *L) {
    const SequencePtr& seq = lua_toseq(L, 1);
    int min = luaL_checkint(L, 2);
    int max = luaL_checkint(L, 3);
    ASSERT_LTE(0, min);
    ASSERT_LTE(min, max);
    ASSERT_LT(max, seq->length());
    int len = max - min + 1;
    const std::string& text = seq->text();
    const char* slice = text.c_str() + min;
    lua_pushlstring(L, slice, len);
    return 1;
}

LUALIB_API int lua_Sequence_tostring(lua_State *L) {
    const SequencePtr& seq = lua_toseq(L, 1);
    std::string repr = seq->tostring();
    lua_pushlstring(L, repr.c_str(), repr.size());
    return 1;
}

LUALIB_API int lua_Sequence_eq(lua_State *L) {
    const SequencePtr& a = lua_toseq(L, 1);
    const SequencePtr& b = lua_toseq(L, 2);
    int eq = ((*a) == (*b));
    lua_pushboolean(L, eq);
    return 1;
}

static const luaL_Reg Sequence_methods[] = {
    {"__gc", lua_Sequence_gc},
    {"type", lua_Sequence_type},
    {"name", lua_Sequence_name},
    {"description", lua_Sequence_description},
    {"genome", lua_Sequence_genome},
    {"chromosome", lua_Sequence_chromosome},
    {"circular", lua_Sequence_circular},
    {"text", lua_Sequence_text},
    {"length", lua_Sequence_length},
    {"sub", lua_Sequence_sub},
    {"__tostring", lua_Sequence_tostring},
    {"__eq", lua_Sequence_eq},
    {NULL, NULL}
};

// -1 is module "model"
static void registerType(lua_State *L,
                         const char* type_name,
                         const char* mt_name,
                         lua_CFunction constructor,
                         const luaL_Reg* methods) {
    luaL_newmetatable(L, mt_name);
    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index"); // mt.__index = mt
    luaL_register(L, NULL, Sequence_methods);
    lua_pushcclosure(L, constructor, 1);
    lua_setfield(L, -2, type_name);
}

extern "C" {
LUALIB_API int luaopen_npge_cmodel(lua_State *L) {
    lua_newtable(L);
    registerType(L, "Sequence", "npge_Sequence",
                 lua_Sequence, Sequence_methods);
    return 1;
}

}
