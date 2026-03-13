/** \file
 * \brief Implementation for pc_tree::PCTree construction methods
 *
 * \author Simon D. Fink <ogdf@niko.fink.bayern>
 *
 * \par License:
 * This file is part of the Open Graph Drawing Framework (OGDF).
 *
 * \par
 * Copyright (C)<br>
 * See README.md in the OGDF root directory for details.
 *
 * \par
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * Version 2 or 3 as published by the Free Software Foundation;
 * see the file LICENSE.txt included in the packaging of this file
 * for details.
 *
 * \par
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * \par
 * You should have received a copy of the GNU General Public
 * License along with this program; if not, see
 * http://www.gnu.org/copyleft/gpl.html
 */

#include "PCNode.h"
#include "PCTree.h"

#include <regex>
#include <stack>

using namespace pc_tree;

PCTree::PCTree(int leafNum, std::vector<PCNode*>* added, PCTreeForest* forest) : PCTree(forest) {
	OGDF_ASSERT(leafNum > 2);
	m_rootNode = newNode(PCNodeType::PNode);
	insertLeaves(leafNum, m_rootNode, added);
}

PCTree::PCTree(const std::string& str, PCTreeNodeArray<std::string>* node_labels, bool keep_ids,
		PCTreeForest* forest)
	: PCTree(forest) {
	if (node_labels) {
		node_labels->init(*this);
	}
	if (keep_ids) { // ensure that all IDs in the string are unique and the nextNodeId is one greater than the max of them
		std::stringstream ss(str);
		std::set<int> ids;
		while (ss.good() && !ss.eof()) {
			if (isdigit(ss.peek())) {
				int idx;
				ss >> idx;
				if (!ids.insert(idx).second) {
					throw std::invalid_argument("Invalid PC-Tree! Illegal re-use of ID "
							+ std::to_string(idx) + " with keep_ids=true at position "
							+ std::to_string(ss.tellg()));
				}
				m_forest->m_nextNodeId = std::max(m_forest->m_nextNodeId, idx + 1);
			} else {
				ss.ignore();
			}
		}
	}

	std::stack<PCNode*> stack;

#define THROW_INVALID                                                                              \
	throw std::invalid_argument("Invalid PC-Tree! Illegal '" + std::string(1, next_char) + "' at " \
			+ "position " + std::to_string(i) + " (parser line " + std::to_string(__LINE__) + ")");

	int i = 0;
	while (i < str.length()) {
		char next_char = str[i];

		if (isspace(next_char)) {
			i++;
			continue;
		}

		PCNode* new_node = nullptr;
		PCNode* parent = stack.empty() ? nullptr : stack.top();
		while (isalnum(next_char) || next_char == '_') {
			std::string next_label;
			while (i < str.length() && (isalnum(next_char) || next_char == '_')) {
				next_label.push_back(next_char);
				i++;
				next_char = str[i];
			}

			int next_id = -1;
			if (keep_ids) {
				next_id = std::stoi(next_label);
				OGDF_ASSERT(next_id >= 0); // '-' would not be alnum
			}

			if (parent == nullptr && getNodeCount() > 0) {
				throw std::invalid_argument("Invalid PC-Tree! Second node '" + next_label + "' at "
						+ "position " + std::to_string(i) + " cannot be at top level (parser line "
						+ std::to_string(__LINE__) + ")");
			}
			new_node = newNode(PCNodeType::Leaf, parent, next_id);
			if (node_labels) {
				(*node_labels)[new_node] = next_label;
			}

			bool delim_seen = false;
			while (i < str.length() && (isspace(next_char) || next_char == ':' || next_char == ',')) {
				if (!isspace(next_char)) {
					if (delim_seen) {
						THROW_INVALID
					}
					delim_seen = true;
					if (next_char == ',') {
						new_node = nullptr; // cannot be used for P-/C-node after comma
					}
				}
				i++;
				next_char = str[i];
			}
		}

		switch (next_char) {
		case '{':
			if (!stack.empty()) {
				THROW_INVALID
			}
			if (!new_node) {
				new_node = newNode(PCNodeType::Leaf);
			}
			if (getLeafCount() != 1 || getNodeCount() != 1) {
				THROW_INVALID
			}
			stack.push(new_node);
			break;
		case '[':
			if (new_node) {
				changeNodeType(new_node, PCNodeType::CNode);
			} else {
				if (parent == nullptr && getNodeCount() > 0) {
					THROW_INVALID // would create a second node at top level
				}
				new_node = newNode(PCNodeType::CNode, parent);
			}
			stack.push(new_node);
			break;
		case '(':
			if (new_node) {
				changeNodeType(new_node, PCNodeType::PNode);
			} else {
				if (parent == nullptr && getNodeCount() > 0) {
					THROW_INVALID // would create a second node at top level
				}
				new_node = newNode(PCNodeType::PNode, parent);
			}
			stack.push(new_node);
			break;
		case ']':
			if (stack.empty() || stack.top()->m_nodeType != PCNodeType::CNode) {
				THROW_INVALID
			}
			stack.pop();
			break;
		case ')':
			if (stack.empty() || stack.top()->m_nodeType != PCNodeType::PNode) {
				THROW_INVALID
			}
			stack.pop();
			break;
		case '}':
			if (stack.empty() || stack.top()->m_nodeType != PCNodeType::Leaf) {
				THROW_INVALID
			}
			stack.pop();
			if (!stack.empty()) {
				THROW_INVALID
			}
			break;
		default:
			THROW_INVALID
		}

		i++;
		next_char = str[i];
		bool delim_seen = false;
		while (i < str.length() && (isspace(next_char) || next_char == ',')) {
			if (!isspace(next_char)) {
				if (delim_seen) {
					THROW_INVALID
				}
				delim_seen = true;
			}
			i++;
			next_char = str[i];
		}
	}
	if (!stack.empty()) {
		throw std::invalid_argument("Invalid PC-Tree! Unexpected end of string");
	}
	OGDF_ASSERT(checkValid());
}

PCTree::PCTree(const PCTree& other, PCTreeNodeArray<PCNode*>& nodeMapping, bool keep_ids,
		PCTreeForest* forest)
	: PCTree(forest) {
	nodeMapping.init(other);
	for (PCNode* other_node : other.allNodes()) {
		PCNode* parent = other_node->getParent();
		int id = -1;
		if (keep_ids) {
			id = other_node->m_id;
			m_forest->m_nextNodeId = std::max(id + 1, m_forest->m_nextNodeId);
		}
		OGDF_ASSERT((parent == nullptr) == (other.m_rootNode == other_node));
		if (parent == nullptr) {
			m_rootNode = nodeMapping[other_node] = newNode(other_node->getNodeType(), nullptr, id);
		} else {
			nodeMapping[other_node] = newNode(other_node->getNodeType(), nodeMapping[parent], id);
		}
	}
	OGDF_ASSERT(nodeMapping[other.m_rootNode] == m_rootNode);
	OGDF_ASSERT(other.getLeafCount() == getLeafCount());
	OGDF_ASSERT(other.getPNodeCount() == getPNodeCount());
	OGDF_ASSERT(other.getCNodeCount() == getCNodeCount());
}

PCTree::~PCTree() {
	OGDF_ASSERT(checkValid());
	// get rid of any degree <= 2 root, including a root leaf and possible degree 2 descendants
	while (m_rootNode != nullptr && m_rootNode->m_childCount <= 2) {
		PCNode* old_root = m_rootNode;
		if (m_rootNode->m_childCount == 2) {
			m_rootNode = old_root->getChild1();
			m_rootNode->detach();
			PCNode* child = old_root->getChild2();
			child->detach();
			m_rootNode->appendChild(child);
		} else if (m_rootNode->m_childCount == 1) {
			m_rootNode = old_root->getOnlyChild();
			m_rootNode->detach();
		} else {
			m_rootNode = nullptr;
		}
		destroyNode(old_root);
	}
	OGDF_ASSERT(checkValid());

	while (!m_leaves.empty()) {
		PCNode* node = m_leaves.back();
		PCNode* parent = node->getParent();
		PCNodeType type = node->getNodeType();
		bool is_root = node == m_rootNode;
		OGDF_ASSERT((parent == nullptr) == is_root);
		if (is_root) {
			OGDF_ASSERT(m_leaves.size() == 1);
			m_rootNode = nullptr;
		}
		node->detach();
		destroyNode(node);
		if (type != PCNodeType::Leaf) {
			m_leaves.pop_back(); // destroyNode(leaf) automatically removes it from the list
		}
		if (parent != nullptr && parent->m_childCount == 0) {
			m_leaves.push_back(parent);
		}
		if (is_root) {
			OGDF_ASSERT(m_leaves.empty());
			OGDF_ASSERT(m_rootNode == nullptr);
		} else {
			OGDF_ASSERT(!m_leaves.empty());
			OGDF_ASSERT(m_rootNode != nullptr);
		}
	}

	if (!m_externalForest) {
		delete m_forest;
	}

	OGDF_ASSERT(m_rootNode == nullptr);
	OGDF_ASSERT(m_pNodeCount == 0);
	OGDF_ASSERT(m_cNodeCount == 0);
}

void PCTree::registerNode(PCNode* node) {
	if (node->m_nodeType == PCNodeType::Leaf) {
		m_leaves.push_back(node);
	} else if (node->m_nodeType == PCNodeType::PNode) {
		m_pNodeCount++;
	} else {
		OGDF_ASSERT(node->m_nodeType == PCNodeType::CNode);
		node->m_nodeListIndex = m_forest->m_parents.makeSet();
		OGDF_ASSERT(m_forest->m_cNodes.size() == node->m_nodeListIndex);
		m_forest->m_cNodes.push_back(node);
		m_cNodeCount++;
	}
}

void PCTree::unregisterNode(PCNode* node) {
	if (node->m_nodeType == PCNodeType::Leaf) {
		m_leaves.erase(node);
	} else if (node->m_nodeType == PCNodeType::PNode) {
		m_pNodeCount--;
	} else {
		OGDF_ASSERT(node->m_nodeType == PCNodeType::CNode);
		OGDF_ASSERT(m_forest->m_cNodes.at(node->m_nodeListIndex) == node);
		m_forest->m_cNodes[node->m_nodeListIndex] = nullptr;
		m_cNodeCount--;
		node->m_nodeListIndex = UNIONFINDINDEX_EMPTY;
	}
}

PCNode* PCTree::newNode(PCNodeType type, PCNode* parent, int id) {
	PCNode* node;
#ifdef OGDF_PCTREE_REUSE_NODES
	if (m_forest->m_reusableNodes) {
		node = m_forest->m_reusableNodes;
		m_forest->m_reusableNodes = m_forest->m_reusableNodes->m_parentPNode;
		node->m_parentPNode = nullptr;
		node->m_timestamp = 0;
		if (id >= 0) {
			node->m_id = id;
			m_forest->m_nextNodeId = std::max(m_forest->m_nextNodeId, id + 1);
		} // else we can't re-use the old ID after clear was called
		else {
			node->m_id = m_forest->m_nextNodeId++;
		}
		node->changeType(type);
	} else
#endif
	{
		if (id < 0) {
			id = m_forest->m_nextNodeId++;
		} else {
			m_forest->m_nextNodeId = std::max(m_forest->m_nextNodeId, id + 1);
		}
		node = new PCNode(m_forest, id, type);
	}
	registerNode(node);
	if (parent != nullptr) {
		parent->appendChild(node);
	} else if (m_rootNode == nullptr) {
		m_rootNode = node;
	}
	m_forest->m_nodeArrayRegistry.keyAdded(node);

	for (auto obs : m_observers) {
		obs->onNodeCreate(node);
	}

	return node;
}

void PCTree::destroyNode(PCNode* const& node) {
	OGDF_ASSERT(node->m_forest == m_forest);
	OGDF_ASSERT(node->isDetached());
	OGDF_ASSERT(node->m_childCount == 0);
	OGDF_ASSERT(node->m_child1 == nullptr);
	OGDF_ASSERT(node->m_child2 == nullptr);
	OGDF_ASSERT(node != m_rootNode);
	unregisterNode(node);
#ifdef OGDF_PCTREE_REUSE_NODES
	node->m_parentPNode = m_forest->m_reusableNodes;
	m_forest->m_reusableNodes = node;
#else
	delete node;
#endif
}

PCNodeType PCTree::changeNodeType(PCNode* node, PCNodeType newType) {
	PCNodeType oldType = node->m_nodeType;
	if (oldType == newType) {
		return oldType;
	}

#ifdef OGDF_DEBUG
	UnionFindIndex oldIndex = node->m_nodeListIndex;
#endif
	unregisterNode(node);
	node->changeType(newType);
	registerNode(node);

	if (oldType == PCNodeType::CNode || newType == PCNodeType::CNode) {
		PCNode* pred = nullptr;
		PCNode* curr = node->m_child1;
#ifdef OGDF_DEBUG
		size_t children = 0;
#endif
		while (curr != nullptr) {
			if (oldType == PCNodeType::CNode) {
				OGDF_ASSERT(curr->m_parentPNode == nullptr);
				OGDF_ASSERT((size_t)m_forest->m_parents.find(curr->m_parentCNodeId) == oldIndex);
			} else {
				OGDF_ASSERT(curr->m_parentPNode == node);
				OGDF_ASSERT(curr->m_parentCNodeId == UNIONFINDINDEX_EMPTY);
			}
			if (newType == PCNodeType::CNode) {
				curr->m_parentPNode = nullptr;
				curr->m_parentCNodeId = node->m_nodeListIndex;
			} else {
				curr->m_parentPNode = node;
				curr->m_parentCNodeId = UNIONFINDINDEX_EMPTY;
			}
#ifdef OGDF_DEBUG
			children++;
#endif
			proceedToNextSibling(pred, curr);
		}
		OGDF_ASSERT(children == node->m_childCount);
		OGDF_ASSERT(pred == node->m_child2);
	}

	return oldType;
}

void PCTree::insertLeaves(int count, PCNode* parent, std::vector<PCNode*>* added) {
	OGDF_ASSERT(parent != nullptr);
	OGDF_ASSERT(parent->m_forest == m_forest);
	if (added) {
		added->reserve(added->size() + count);
	}
	for (int i = 0; i < count; i++) {
		PCNode* leaf = newNode(PCNodeType::Leaf);
		parent->appendChild(leaf);
		if (added) {
			added->push_back(leaf);
		}
	}
}

void PCTree::replaceLeaf(int leafCount, PCNode* leaf, std::vector<PCNode*>* added) {
	OGDF_ASSERT(leaf && leaf->m_forest == m_forest);
	OGDF_ASSERT(leaf->isLeaf());
	OGDF_ASSERT(leafCount > 1);
	if (getLeafCount() <= 2) {
		changeNodeType(leaf->getParent(), PCNodeType::PNode);
		insertLeaves(leafCount, leaf->getParent(), added);
		leaf->detach();
		destroyNode(leaf);
	} else {
		changeNodeType(leaf, PCNodeType::PNode);
		insertLeaves(leafCount, leaf, added);
	}
}

void PCTree::destroyLeaf(PCNode* leaf) {
	OGDF_ASSERT(leaf->getNodeType() == PCNodeType::Leaf);
	OGDF_ASSERT(leaf != m_rootNode);

	PCNode* parent = leaf->getParent();
	leaf->detach();
	destroyNode(leaf);

	/* assume the PC-tree is valid, so a childCount of 0 is impossible, since every inner node must
     * have at least 2 children (except for child of root node) */

	if (parent->m_childCount == 1) {
		if (m_rootNode->getNodeType() == PCNodeType::Leaf) {
			if (parent->getChild1()->getNodeType() != PCNodeType::Leaf
					|| m_rootNode->getChild1() != parent) {
				PCNode* child = parent->getChild1();
				child->detach();
				parent->replaceWith(child);
				destroyNode(parent);
			}
		} else {
			PCNode* child = parent->getChild1();
			if (parent != m_rootNode) {
				child->detach();
				parent->replaceWith(child);
				destroyNode(parent);
			} else if (child->getNodeType() != PCNodeType::Leaf) {
				PCNode* root = m_rootNode;
				root->detach();
				root->m_childCount = 0;
				root->m_child1 = root->m_child2 = nullptr;
				child->m_parentCNodeId = UNIONFINDINDEX_EMPTY;
				child->m_parentPNode = nullptr;
				setRoot(child);
				destroyNode(root);
			}
		}
	}
}

void PCTree::insertTree(PCNode* at, PCTree* inserted) {
	OGDF_HEAVY_ASSERT(checkValid());
	OGDF_HEAVY_ASSERT(inserted->checkValid());
	OGDF_ASSERT(at->isValidNode(getForest()));
	OGDF_ASSERT(inserted->getForest() == getForest());

	m_observers.splice(m_observers.end(), inserted->m_observers);
	m_leaves.splice(m_leaves.end(), inserted->m_leaves);
	m_pNodeCount += inserted->m_pNodeCount;
	m_cNodeCount += inserted->m_cNodeCount;

	PCNode* root = inserted->m_rootNode;
	inserted->m_rootNode = nullptr;
	inserted->m_pNodeCount = inserted->m_cNodeCount = 0;
	delete inserted;
	inserted = nullptr;

	while (root->getChildCount() == 1) {
		PCNode* new_root = root->getOnlyChild();
		new_root->detach();
		destroyNode(root);
		root = new_root;
	}

	if (at->isLeaf()) {
		OGDF_ASSERT(!at->isDetached());
		PCNode* parent = at->getParent();
		if (!root->isLeaf() && parent->getDegree() == 2) {
			// remove degree 2 node in a tree with two leaves if it got further leaves
			at->detach();
			if (parent->isDetached()) {
				// use the inserted root as root
				m_rootNode = root;
				PCNode* child = parent->getOnlyChild();
				child->detach();
				m_rootNode->appendChild(child);
			} else {
				// we can use the other leaf as root
				parent->replaceWith(root);
			}
			destroyNode(parent);
		} else {
			at->replaceWith(root);
		}
		destroyNode(at);
	} else {
		at->appendChild(root);
	}
	OGDF_ASSERT(checkValid());
}
