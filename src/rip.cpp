#include "rip.h"
// #include "set_ops.h"  // For set_any
// #include <unordered_map>

// using VS = std::vector<std::string>;

//' Maximum Cardinality Search
//' 
//' @param adj A named adjacency list of a decomposable graph
//' @param check Boolean: check if adj is decomposable
//' @details If adj is not the adjacency list of a decomposable graph an error is raised
//' @return A list with a perfect numbering of the nodes and a perfect sequence of sets
//' @examples
//' x <- list(a = c("b", "d"), b = c("a", "c", "d"), c = c("b", "d"), d = c("a", "c", "b"))
//' mcs(x)
//' @export
// [[Rcpp::export]]
Rcpp::List mcs(Rcpp::List & adj, bool check = true) {
  VS  nodes = adj.names();
  int N = nodes.size();
  // Raise a WARNING if adj is the empty graph
  // if( !( N - 1 )) return VS(nodes[0]);
  std::unordered_map<std::string, int> labels; // = {};
  for( int i = 0; i < N; i++ ) {
    labels.emplace(nodes[i], 0);
  }
  VVS ps(N);
  decltype(nodes) remaining_nodes = nodes;
  decltype(nodes) used_nodes(N, "");
  auto v = nodes[0];
  used_nodes[0] = v;
  ps[0] = {v};
  remaining_nodes.erase(remaining_nodes.begin()+0);
  for( int i = 1; i < N; i++ ) {
    auto ne_i = Rcpp::as<VS>(adj[v]);
    // Increment neighbor nodes with a one
    for (auto it = ne_i.begin(); it != ne_i.end(); ++it) {
      auto ne_ = labels.find(*it);
      ne_->second++;
    }
    std::string max_v;
    int max_val = -1;
    VS::iterator max_it;
    for (auto it = remaining_nodes.begin(); it != remaining_nodes.end(); ++it) {
      auto rn = labels.find(*it);
      int max_candidate = rn->second;
      if( max_candidate > max_val ) {
	max_v   = *it;
	max_val = max_candidate;
	max_it  = it;
      }
    }
    v = max_v;
    used_nodes[i] = v;
    remaining_nodes.erase(max_it);
    auto ne_v = Rcpp::as<VS>(adj[v]);
    ne_v.push_back(v); // The closure of v
    VS anc = VS(used_nodes.begin(), used_nodes.begin() + i + 1);
    VS B_i = set_intersect(ne_v, anc);
    if (check) {
      int card_i = B_i.size();
      if ( i > 1 && card_i > 2) {
      // ----------------------------------------------------------------------------------------------
      // Test for decomposability for step i. See Lauritzen for details
      // 1. cl(v_i) \cap {v_1, .., v_{i-1}} needs to be complete
      // 2. The check is always positive for i \in {1,2}
      // 3. It is not necessary to check for ||ne_v|| < 3 since these are always complete in the graph
      // ----------------------------------------------------------------------------------------------
	for (int j = 0; j < card_i; j++) {
	  for (int k = j + 1; k < card_i; k++) { // k < card_i - 1
	    auto adj_k = Rcpp::as<VS>(adj[B_i[k]]);
	    if ( !set_in(B_i[j], adj_k) ) Rcpp::stop("The corresponding graph of <adj> is not decomposable");
    	}
      }
     }
    }
    ps[i] = B_i;
  }
  Rcpp::List ps_out = Rcpp::wrap(ps);
  ps_out.names() = used_nodes;
  return Rcpp::List::create(Rcpp::_["po"] = used_nodes , Rcpp::_["ps"] = ps_out);
} 

// [[Rcpp::export]]
VVS perfect_cliques(VVS & x) {
  // In: 
  // x: a perfect sequence of sets (denoted B_i in Lauritzen)
  //
  // Out:
  // y: a perfect sequence of the cliques
  int n = x.size();
  VVS pc;
  for(int i = 0; i < n; i++) {
    std::vector<bool> v; // Dummy var to loop over
    for(int j = 0; j < n; j++) {
      if( j != i ) {
	v.push_back(set_issubeq(x[i], x[j]));
      }
    }
    if( !set_any(v) ) {
      pc.push_back(x[i]); 
    }
  }
  return pc;
}

// [[Rcpp::export]]
Rcpp::List perfect_separators(VVS & x) {
  // x: Cliques (with RIP ordering)
  // S_j := H_{j-1} \cap C_j, H_{j-1} := \cap_k C_{k-1}, k = 1, 2, ..., j-1
  int n = x.size();
  Rcpp::List ps(n);       // All elements initialized to NULL
  if( n == 1 ) return ps; // List::create(_[""] = R_NilValue);
  for (int i = 1; i < n; i++) {
    VS Hi_1;
    for (int j = 0; j < i; j++) {
      Hi_1.insert(Hi_1.end(), x[j].begin(), x[j].end());
    }
    ps[i] = set_intersect(x[i], Hi_1);
  }
  return ps;
}

// [[Rcpp::export]]
Rcpp::List parents(VS po, Rcpp::List ps) {
  // ps: perfect sequence from mcs_marked
  int npo = po.size();
  for (int i = 0; i < npo; i++) {
    std::string poi = po[i];
    auto psi    = Rcpp::as<VS>(ps[i]);
    auto psi_it = std::find(psi.begin(), psi.end(), poi);
    if (psi_it != psi.end()) {
      psi.erase(psi_it);
      ps[i] = psi;
    }
  }
  return ps;
}

//' Running Intersection Property
//' @description Given a decomposable graph, this functions finds a perfect numbering on the vertices using maximum cardinality search, and hereafter returns a list with two elements: "C" - the cliques and "S" - the separators.
//'
//' @param adj A named adjacency list of a decomposable graph
//' @param check Boolean: check if adj is decomposable
//' @seealso \code{\link{mcs}}, \code{\link{is_decomposable}}
//' @description Given a decomposable graph, this functions finds a perfect numbering on the vertices
//' using maximum cardinality search, and hereafter returns a list with three elements:
//' "C"   - the cliques,
//' "S"   - the separators and
//' "PA" - the parents (all vertices adjacent to the node i that are numbered before node i)
//' @examples
//' x <- list(a = c("b", "d"), b = c("a", "c", "d"), c = c("b", "d"), d = c("a", "c", "b"))
//' y <- rip(x)
//' # Cliques:
//' y$C
//' # Separators:
//' y$S
//' # Parents:
//' y$PA
//' @export
// [[Rcpp::export]]
Rcpp::List rip(Rcpp::List & adj, bool check = true) {
  auto  z       = mcs(adj, check);
  VVS pseq      = z["ps"];
  VVS pc        = perfect_cliques(pseq);
  Rcpp::List ps = perfect_separators(pc);
  Rcpp::List pa = parents(z["po"], z["ps"]);
  return Rcpp::List::create(Rcpp::_["C"] = pc , Rcpp::_["S"] = ps, Rcpp::_["PA"] = pa);
}


/*----------------------------------------------------------*
 *                     Marked Graphs                                   
 * ---------------------------------------------------------*/
void insert_star(Rcpp::List & adj,  VS & disc_vars) {
  adj.insert(adj.begin(), disc_vars);
  VS new_names = adj.names();
  new_names[0] = "_*star*_";
  adj.names() = new_names;
  for (auto & e : disc_vars) {
    VS adje = adj[e];
    adje.insert(adje.begin(), "_*star*_");
    adj[e] = adje;
  }
  disc_vars.insert(disc_vars.begin(), "_*star*_");
}
void remove_star(VVS & ps, VS & used_nodes, int nd) {
  used_nodes.erase(used_nodes.begin()); // begin is/should be the star
  ps.erase(ps.begin());
  for (int i = 0; i < nd; i++) {
    // The first |disc_vars| elements _ARE/SHOULD_ be the discrete ones!
    auto & psi = ps[i];
    auto it = std::find(psi.begin(), psi.end(), "_*star*_");
    if (it != psi.end()) {
      psi.erase(it); 
    }
  }
}

//' Maximum Cardinality Search in Marked Graphs
//' 
//' @param adj A named adjacency list of a decomposable marked graph
//' @param disc_vars Character vector of discrete variables
//' @param check Boolean: check if adj is decomposable
//' @details If adj is not the adjacency list of a decomposable marked graph an error is raised
//' @return A list with a perfect numbering of the nodes and a perfect sequence of sets
//' @examples
//' x <- list(a = c("b", "d"), b = c("a", "c", "d"), c = c("b", "d"), d = c("a", "c", "b"))
//' mcs_marked(x, c("a"))
//' @export
// [[Rcpp::export]]
Rcpp::List mcs_marked(Rcpp::List & adj, VS disc_vars, bool check = true) {
  if (check) {
    insert_star(adj, disc_vars);
  }
  VS nodes = adj.names();
  int N = nodes.size();
  // Raise a WARNING if adj is the empty graph
  // if( !( N - 1 )) return VS(nodes[0]);
  std::unordered_map<std::string, int> labels; // = {};
  for( int i = 0; i < N; i++ ) {
    labels.emplace(nodes[i], 0);
  }
  VVS ps(N);
  decltype(nodes) remaining_nodes = nodes;
  decltype(nodes) used_nodes(N, "");
  auto v = disc_vars[0]; // Current vertex
  used_nodes[0] = v;
  ps[0] = {v};
  auto rm_v = std::find(remaining_nodes.begin(), remaining_nodes.end(), v);
  remaining_nodes.erase(rm_v);
  for( int i = 1; i < N; i++ ) {
    auto ne_i = Rcpp::as<VS>(adj[v]);
    // Increment neighbor nodes with a one
    for (auto it = ne_i.begin(); it != ne_i.end(); ++it) {
      auto ne_ = labels.find(*it);
      ne_->second++;
    }
    std::string max_v;
    int max_val = -1;
    VS::iterator max_it;
    VS remaining_discrete_nodes = set_intersect(disc_vars, remaining_nodes);
    VS remaining_container;
    if (!remaining_discrete_nodes.empty()) {
      remaining_container = remaining_discrete_nodes;
    } else {
      remaining_container = remaining_nodes;
    }
    for (auto it = remaining_container.begin(); it != remaining_container.end(); ++it) {
      auto rn = labels.find(*it);
      int max_candidate = rn->second;
      if( max_candidate > max_val ) {
	max_v   = *it;
	max_val = max_candidate;
	max_it  = it;
      }
    }
    v = max_v;
    used_nodes[i] = v;
    auto rm_v = std::find(remaining_nodes.begin(), remaining_nodes.end(), v);
    remaining_nodes.erase(rm_v);
    auto ne_v = Rcpp::as<VS>(adj[v]);
    ne_v.push_back(v); // The closure of v.
    VS anc = VS(used_nodes.begin(), used_nodes.begin() + i + 1);
    VS B_i = set_intersect(ne_v, anc); // B_i sets in Laurtizen 96.
    if (check) {
      // ----------------------------------------------------------------------------------------------
      // Test for decomposability for step i. See Lauritzen for details
      // 1. cl(v_i) \cap {v_1, .., v_{i}} needs to be complete
      // 2. The check is always positive for i \in {1,2}
      // 3. It is not necessary to check for ||ne_v|| < 3 since these are always complete in the graph
      // ----------------------------------------------------------------------------------------------      
      int card_i = B_i.size();
      if (i > 1 && card_i > 2) {
	for (int j = 0; j < card_i; j++) {
	  for (int k = j + 1; k < card_i; k++) {
	    auto adj_k = Rcpp::as<VS>(adj[B_i[k]]);
	    if (!set_in(B_i[j], adj_k)) Rcpp::stop("The corresponding graph of <adj> is not decomposable");
    	  }
        }
      }
    }
    ps[i] = B_i;
  }
  if (check) { //remove star from "po" and "ps"
    remove_star(ps, used_nodes, disc_vars.size() - 1); // disc_vars has star as first element
  }
  Rcpp::List ps_out = Rcpp::wrap(ps);
  ps_out.names() = used_nodes;
  return Rcpp::List::create(Rcpp::_["po"] = used_nodes , Rcpp::_["ps"] = ps_out);
}


//' Running Intersection Property for Marked Graphs
//' @description Given a decomposable marked graph, this functions finds a perfect numbering on the vertices
//' using maximum cardinality search, and hereafter returns a list with three elements:
//' "C"   - the cliques,
//' "S"   - the separators and
//' "PA" - the parents (all vertices adjacent to the node i that are numbered before node i)
//'
//' @param adj A named adjacency list of a decomposable marked graph
//' @param disc_vars Character vector of discrete variables in adj
//' @param check Boolean: check if adj is decomposable
//' @seealso \code{\link{mcs_marked}}, \code{\link{is_decomposable}} 
//' @examples
//' x <- list(a = c("b", "d"), b = c("a", "c", "d"), c = c("b", "d"), d = c("a", "c", "b"))
//' y <- rip_marked(x, c("a"))
//' # Cliques:
//' y$C
//' # Separators:
//' y$S
//' # Parents:
//' y$PA
//' @export
// [[Rcpp::export]]
Rcpp::List rip_marked(Rcpp::List & adj, VS disc_vars, bool check = true) {
  auto z         = mcs_marked(adj, disc_vars, check);
  VVS zps        = z["ps"];
  VVS pc         = perfect_cliques(zps);
  Rcpp::List ps  = perfect_separators(pc);
  Rcpp::List pa  = parents(z["po"], z["ps"]);
  return Rcpp::List::create(Rcpp::_["C"] = pc , Rcpp::_["S"] = ps, Rcpp::_["PA"] = pa);
}
