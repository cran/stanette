#ifndef STAN_LANG_AST_FUN_GENERATE_EXPRESSION_HPP
#define STAN_LANG_AST_FUN_GENERATE_EXPRESSION_HPP

#include <stan/lang/ast.hpp>
#include <ostream>

namespace stan {
namespace lang {

struct expression;

/**
 * Write the code generated by the specified expression to the
 * specified output stream, putting it in a user-readable format
 * if the user-facing flag is true.  This is just the header for a
 * forward declaration defined in the generator.
 *
 * @param[in] e expression to write
 * @param[in] user_facing true if expression should be written so
 * that a user can understand it
 * @param[in, out] o stream to which expression is written
 */
void generate_expression(const expression& e, bool user_facing,
                         std::ostream& o);
}  // namespace lang
}  // namespace stan
#endif
