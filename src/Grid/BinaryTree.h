/* 
 * File:   BinaryTree.h
 * Author: Orhan Shibliyev
 *
 * Created on June 13, 2014, 5:51 AM
 */

#ifndef BINARYTREE_H
#define	BINARYTREE_H

#include <map>

struct bnode
{
	int key_value;
	bnode *left;
	bnode *right;
};

class btree
{
public:

	std::vector <unsigned int> a;

	bnode *root;
	btree ()
	{
		root = NULL;
	}
	~btree ()
	{
		destroy_tree ();
	}

	void insert (int key)
	{
		if (root != NULL)
		{
			insert (key, root);
		}
		else
		{
			root = new bnode;
			root -> key_value = key;
			root -> left = NULL;
			root -> right = NULL;
			a.push_back (key);
		}
	}
	bool search (int key)
	{
		return search (key, root);
	}
	void destroy_tree ()
	{
		destroy_tree (root);
	}

private:

	void destroy_tree (bnode *&leaf)
	{
		if (leaf != NULL)
		{
			destroy_tree (leaf -> left);
			destroy_tree (leaf -> right);
			delete leaf;
			leaf = NULL;
		}
	}

	void insert (int key, bnode *leaf)
	{
		if (key < leaf -> key_value)
		{
			if (leaf -> left != NULL)
			{
				insert (key, leaf -> left);
			}
			else
			{
				leaf -> left = new bnode;
				leaf -> left -> key_value = key;
				leaf -> left -> left = NULL;
				leaf -> left -> right = NULL;
				a.push_back (key);
			}
		}
		else if (key >= leaf -> key_value)
		{
			if (leaf -> right != NULL)
			{
				insert (key, leaf -> right);
			}
			else
			{
				leaf -> right = new bnode;
				leaf -> right -> key_value = key;
				leaf -> right -> left = NULL;
				leaf -> right -> right = NULL;
				a.push_back (key);
			}
		}
	}

	bool search (int key, bnode *leaf)
	{
		if (leaf != NULL)
		{
			if (key == leaf -> key_value)
			{
				return true;
			}
			if (key < leaf -> key_value)
			{
				return search (key, leaf -> left);
			}
			else
			{
				return search (key, leaf -> right);
			}
		}
		else
		{
			return NULL;
		}
	}
};

struct faceNumCheck
{
	std::map <int, bool> nmap;
};

#endif	/* BINARYTREE_H */

