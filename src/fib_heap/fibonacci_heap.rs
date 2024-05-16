use std::collections::VecDeque;

#[derive(Debug, Clone)]
pub struct InternalTree<T: std::cmp::Ord, K> {
    // number of direct children of the current node
    degree: usize,

    // data stored in the current node
    payload: (T, K),
    // children of the current node
    children_list: VecDeque<InternalTree<T, K>>,

    // indicates wether current node is a min heap-ordered tree or not
    min: bool,
}

impl<T: std::cmp::Ord, K> InternalTree<T, K> {
    // initializes an internal tree which is min or max heap ordered tree based on min parameter
    fn init(payload: (T, K), min: bool) -> InternalTree<T, K> {
        InternalTree {
            degree: 0,
            payload,
            children_list: VecDeque::new(),
            min,
        }
    }

    // returns true if tree1.payload <= tree2.payload
    #[inline]
    fn is_smaller_or_equal(
        internal_tree_1: &InternalTree<T, K>,
        internal_tree_2: &InternalTree<T, K>,
    ) -> bool {
        let payload1 = internal_tree_1.peek_payload();
        let payload2 = internal_tree_2.peek_payload();
        payload1.0 <= payload2.0
    }
    // returns true if tree1.payload >= tree2.payload
    #[inline]
    fn is_greater_or_equal(
        internal_tree_1: &InternalTree<T, K>,
        internal_tree_2: &InternalTree<T, K>,
    ) -> bool {
        let payload1 = internal_tree_1.peek_payload();
        let payload2 = internal_tree_2.peek_payload();
        payload1.0 >= payload2.0
    }

    // returns true if tree1 has higher priority than tree2
    // it means:
    // if trees are min heap-ordered higher priority means smaller values
    // if trees are max heap-ordered higher priority means larger values
    fn has_higher_priority(
        internal_tree_1: &InternalTree<T, K>,
        internal_tree_2: &InternalTree<T, K>,
        is_min: bool,
    ) -> bool {
        if is_min {
            InternalTree::is_smaller_or_equal(&internal_tree_1, &internal_tree_2)
        } else {
            InternalTree::is_greater_or_equal(&internal_tree_1, &internal_tree_2)
        }
    }

    // merges two heap-ordered trees and returns the merged tree
    fn merge(
        mut internal_tree_1: InternalTree<T, K>,
        mut internal_tree_2: InternalTree<T, K>,
    ) -> InternalTree<T, K> {
        // make sure both tree are of the same kind
        if internal_tree_1.is_min() != internal_tree_2.is_min() {
            panic!("Both internal trees must be of same type. Both min or both max")
        }

        let trees_are_min = internal_tree_1.is_min();

        // tree with lower priority must be child of the tree with higher priority
        if InternalTree::has_higher_priority(&internal_tree_1, &internal_tree_2, trees_are_min) {
            internal_tree_1.add_child(internal_tree_2);

            internal_tree_1
        } else {
            internal_tree_2.add_child(internal_tree_1);

            internal_tree_2
        }
    }

    // add another internal tree as a child
    fn add_child(&mut self, internal_tree: InternalTree<T, K>) {
        self.children_list.push_back(internal_tree);
        self.degree += 1;
    }

    // returns degree of current tree
    fn degree(&self) -> usize {
        self.degree
    }

    // returns a reference to root payload of current tree
    fn peek_payload(&self) -> &(T, K) {
        &self.payload
    }

    // returns a reference to list of children of the current node
    fn children_list(&self) -> &VecDeque<InternalTree<T, K>> {
        &self.children_list
    }

    // returns true if tree is initialized as a min heap-ordered tree
    fn is_min(&self) -> bool {
        self.min
    }
}

// ------------- Fibonacci Heap -------------
/// A Fibonacci heap is a data structure for priority queue operations.
/// It has a better amortized running time than binary heap and binomial heap.
///
/// # Examples
/// ```
/// use rudac::heap::FibonacciHeap;
///
/// // initialize a fibonacci heap
/// let mut fibonacci_heap: FibonacciHeap<usize> = FibonacciHeap::init_min();
///
/// // push items into heap
/// fibonacci_heap.push(0);
/// fibonacci_heap.push(1);
/// fibonacci_heap.push(3);
///
/// // heap will have the shape:
/// //  min
/// //   |
/// //   0 <-> 1 <-> 2
/// assert_eq!(
///     FibonacciHeap::preorder(&fibonacci_heap),
///     String::from("Priority: 0\nTree 1: 1\nTree 2: 3\n")
/// )
/// ```
#[derive(Debug, Clone)]
pub struct FibonacciHeap<T: std::cmp::Ord, K> {
    // doubly linked list of internal trees
    children_list: VecDeque<InternalTree<T, K>>,

    // total number of items in the heap
    size: usize,

    // pointer to root containing the highest priority
    priority_pointer: Option<InternalTree<T, K>>,

    // indicates wether current heap is initialized as a min heap or not
    min: bool,
}

impl<T: std::cmp::Ord, K> FibonacciHeap<T, K>
where
    T: std::clone::Clone,
    K: std::clone::Clone,
    (T, K): std::cmp::PartialEq,
{
    // initializes a fibonacci heap
    fn init(min: bool) -> FibonacciHeap<T, K> {
        FibonacciHeap {
            children_list: VecDeque::new(),
            size: 0,
            priority_pointer: None,
            min,
        }
    }

    /// Initializes a min heap with the specified `payload`
    ///
    /// # Examples
    /// ```
    /// use rudac::heap::FibonacciHeap;
    ///
    /// let fibonacci_heap: FibonacciHeap<usize> = FibonacciHeap::init_min();
    ///
    /// assert_eq!(fibonacci_heap.is_min(), true);
    /// ```
    pub fn init_min() -> FibonacciHeap<T, K> {
        FibonacciHeap::init(true)
    }

    /// Initializes a max heap with the specified `payload`
    ///
    /// # Examples
    /// ```
    /// use rudac::heap::FibonacciHeap;
    ///
    /// let fibonacci_heap: FibonacciHeap<usize> = FibonacciHeap::init_max();
    ///
    /// assert_eq!(fibonacci_heap.is_max(), true);
    /// ```
    pub fn init_max() -> FibonacciHeap<T, K> {
        FibonacciHeap::init(false)
    }

    /// Pushes specified `payload` into heap
    ///
    /// # Arguments:
    /// * `payload`: data to be pushed into heap
    ///
    /// # Examples
    /// ```
    /// use rudac::heap::FibonacciHeap;
    ///
    /// let mut fibonacci_heap: FibonacciHeap<usize> = FibonacciHeap::init_min();
    ///
    /// // push items into heap
    /// fibonacci_heap.push(0);
    /// fibonacci_heap.push(1);
    /// fibonacci_heap.push(3);
    ///
    /// // heap will have the shape:
    /// //  min
    /// //   |
    /// //   0 <-> 1 <-> 2
    /// assert_eq!(
    ///     FibonacciHeap::preorder(&fibonacci_heap),
    ///     String::from("Priority: 0\nTree 1: 1\nTree 2: 3\n")
    /// )
    /// ```
    pub fn push(&mut self, payload: (T, K)) {
        // create a compatible root with current heap, containing the payload
        let new_node = InternalTree::init(payload, self.is_min());

        let heap_is_min = self.is_min();

        // if there is no priority node, assign the newly created node as priority node
        if self.priority_pointer.is_none() {
            self.priority_pointer = Some(new_node);
        } else {
            if InternalTree::has_higher_priority(
                // if new node has higher priority, it must become priority node
                &new_node,
                &self.priority_pointer.as_ref().unwrap(),
                heap_is_min,
            ) {
                // swap new node and priority node
                let temp = self.priority_pointer.take().unwrap();
                self.priority_pointer = Some(new_node);
                self.children_list.push_back(temp);
            } else {
                // if new node has lower priority, just add it to children list of the heap
                self.children_list.push_back(new_node);
            }
        }

        // account for newly added node
        self.size += 1;
    }

    /// Merges two fibonacci heaps and returns the merged fibonacci heap
    ///
    /// # Arguments:
    /// * `fibonacci_heap_1`: first fibonacci heap
    /// * `fibonacci_heap_2`: second fibonacci heap
    ///
    /// # Panics:
    /// * panics if two fibonacci heaps are not the same kind(ex. one is min heap and the other is max heap)
    ///
    /// # Examples
    /// ```
    /// use rudac::heap::FibonacciHeap;
    ///
    /// let mut fibonacci_heap_1: FibonacciHeap<usize> = FibonacciHeap::init_min();
    /// fibonacci_heap_1.push(0);
    /// fibonacci_heap_1.push(2);
    ///
    /// let mut fibonacci_heap_2: FibonacciHeap<usize> = FibonacciHeap::init_min();
    /// fibonacci_heap_2.push(1);
    /// fibonacci_heap_2.push(3);
    ///
    /// let merged_heap = FibonacciHeap::merge(fibonacci_heap_2, fibonacci_heap_1);
    ///
    /// assert_eq!(
    ///     FibonacciHeap::preorder(&merged_heap),
    ///     String::from("Priority: 0\nTree 1: 3\nTree 2: 2\nTree 3: 1\n")
    /// );
    /// ```
    pub fn merge(
        mut fibonacci_heap_1: FibonacciHeap<T, K>,
        mut fibonacci_heap_2: FibonacciHeap<T, K>,
    ) -> FibonacciHeap<T, K> {
        // if one heap is min and the other is max, panic!. merge is not possible
        if fibonacci_heap_1.is_min() != fibonacci_heap_2.is_min() {
            panic!("Two heaps must be of same type in order for merge to be possible")
        }

        // if either heaps are empty, return the other one as result
        if fibonacci_heap_1.is_empty() {
            return fibonacci_heap_2;
        } else if fibonacci_heap_2.is_empty() {
            return fibonacci_heap_1;
        }

        // concatenate children list of heap1 and heap2
        fibonacci_heap_1
            .children_list
            .append(&mut fibonacci_heap_2.children_list);

        let heap_is_min = fibonacci_heap_1.is_min();

        // update priority node in merged heap
        // if priority node in heap2 has higher priority than priority node in heap1, priority node in heap2 must become the new priority node of merged heap
        if InternalTree::has_higher_priority(
            &fibonacci_heap_2.priority_pointer.as_ref().unwrap(),
            &fibonacci_heap_1.priority_pointer.as_ref().unwrap(),
            heap_is_min,
        ) {
            // swap priority nodes of heap2 and heap1
            let temp = fibonacci_heap_1.priority_pointer.take().unwrap();
            fibonacci_heap_1.priority_pointer = fibonacci_heap_2.priority_pointer.take();
            fibonacci_heap_1.children_list.push_back(temp);
        } else {
            // if priority node of heap2 has lower priority then just add it to children list of heap1
            fibonacci_heap_1
                .children_list
                .push_back(fibonacci_heap_2.priority_pointer.unwrap());
        }

        // calculate size of merged heap
        fibonacci_heap_1.size += fibonacci_heap_2.size;

        // return merged heap
        fibonacci_heap_1
    }

    /// Pops and returns item with highest priority. Returns `None` if heap is empty. After pop, heap will be consolidated
    ///
    /// # Examples
    /// ```
    /// use rudac::heap::FibonacciHeap;
    ///
    /// let mut fibonacci_heap: FibonacciHeap<usize> = FibonacciHeap::init_min();
    /// fibonacci_heap.push(2);
    /// fibonacci_heap.push(3);
    /// fibonacci_heap.push(0);
    /// fibonacci_heap.push(1);
    ///
    /// // before pop
    /// assert_eq!(
    ///     FibonacciHeap::preorder(&fibonacci_heap),
    ///     String::from("Priority: 0\nTree 1: 3\nTree 2: 2\nTree 3: 1\n")    
    /// );
    ///
    /// assert_eq!(fibonacci_heap.pop(), Some(0));
    ///
    /// // heap trees are consolidated
    /// assert_eq!(
    ///     FibonacciHeap::preorder(&fibonacci_heap),
    ///     String::from("Priority: 1\nTree 1: 2 3\n")
    /// );
    /// ```
    pub fn pop(&mut self) -> Option<(T, K)> {
        if self.is_empty() {
            return None;
        }

        // extract node with highest priority from heap
        let priority_node = self.priority_pointer.take().unwrap();

        // account for deleted node
        self.size -= 1;

        // iterate over children of removed node and add them to children list of heap
        self.children_list.extend(priority_node.children_list);

        // extract payload of priority node
        let payload = priority_node.payload;

        // if there is nodes in heap, consolidate them
        if !self.is_empty() {
            // a temp priority node just for consolidate method to work
            self.priority_pointer = self.children_list.pop_front();

            self.consolidate();
        }

        // return payload with highest priority
        Some(payload)
    }

    // this method consolidate trees in fibonacci heap
    // until each tree in children list of the heap has a unique degree
    // ex after consolidate there can not be two trees with degree 0 like: 0 <-> 1
    //
    // after consolidate total number of trees in heap is gonna be at most log(n)
    fn consolidate(&mut self) {
        // there is nothing to consolidate
        if self.is_empty() {
            return;
        }
        // use a helper vector for consolidating
        // vector keeps track of degree of present trees
        // therefore we can make sure each degree is associated with a unique tree
        // array size will be log(heap size) with base 1.61803
        let array_size = ((self.size as f32).log(1.61803_f32) + 1.0) as usize;

        // helper vector for tracking current degrees present in consolidating process
        let mut a: Vec<Option<InternalTree<T, K>>> = Vec::with_capacity(array_size);

        // initialize consolidate array
        a.resize_with(array_size, || None);

        // add priority node to children list
        // because we have to iterate over all nodes
        self.children_list
            .push_front(self.priority_pointer.take().unwrap());

        // iterate over children and merge trees with same degrees
        for mut x in self.children_list.drain(..) {
            let mut d = x.degree(); // degree of current internal tree
            while a[d].is_some() {
                // iterate over consolidate array to find the place for x
                let y = a[d].take().unwrap(); // if there exists a tree with degree of x like y
                x = InternalTree::merge(x, y); // merge x and y and store merged tree in x
                d += 1; // degree of x is now d + 1 because it has y as its child
            }
            a[d] = Some(x); // finally when a degree is free(a[d]), means degree of x is unique. store it in consolidate array
        }

        // update priority pointer and children list
        let heap_is_min = self.is_min();

        // after consolidate, "a" has all the nodes in the heap
        // we have to find minimum between these nodes and add rest of them to children list of heap
        // so iterate over consolidate array

        let mut nodes = a.into_iter().filter_map(|x| x);
        let mut priority_pointer = nodes.next().unwrap();

        for mut node in nodes {
            if InternalTree::has_higher_priority(&node, &priority_pointer, heap_is_min) {
                // current tree in a has higher priority than latest found priority node, swap them
                std::mem::swap(&mut priority_pointer, &mut node);
            }
            self.children_list.push_back(node);
        }

        self.priority_pointer = Some(priority_pointer);
    }

    pub fn decrease_key(&mut self, old_payload: (T, K), new_payload: (T, K)) -> bool {
        let mut node_index = None;

        // Find the index of the node with the old key
        for (index, internal_tree) in self.children_list.iter().enumerate() {
            if internal_tree.peek_payload() == &old_payload {
                node_index = Some(index);
                break;
            }
        }

        // If the node with the old key is found, decrease its key
        if let Some(index) = node_index {
            let mut temp_heap = FibonacciHeap::init(self.is_min()); // Create a temporary heap

            // Split the heap into two parts: nodes before the node with the old key and after it
            let mut temp_children_list = self.children_list.split_off(index + 1);
            std::mem::swap(&mut self.children_list, &mut temp_children_list); // Swap back the original children list

            let mut node_to_decrease = temp_children_list.pop_front().unwrap(); // Remove the node with the old key
            node_to_decrease.payload = new_payload; // Decrease its key

            // Merge the node with the decreased key back into the heap
            temp_heap.children_list = temp_children_list;
            temp_heap.priority_pointer = Some(node_to_decrease);
            *self = FibonacciHeap::merge(self.clone(), temp_heap);

            true // Key successfully decreased
        } else {
            false // Node with the old key not found
        }
    }

    /// Returns a reference to item with highest priority
    ///
    /// # Examples
    /// ```
    /// use rudac::heap::FibonacciHeap;
    ///
    /// let mut fibonacci_heap: FibonacciHeap<usize> = FibonacciHeap::init_min();
    ///
    /// fibonacci_heap.push(0);
    ///
    /// assert_eq!(fibonacci_heap.peek(), Some(&0));
    /// ```
    pub fn peek(&self) -> Option<&(T, K)> {
        if self.is_empty() {
            return None;
        }

        let payload = self.priority_pointer.as_ref().unwrap().peek_payload();
        Some(payload)
    }

    /// Clears the heap and resets internal flags
    ///
    /// # Examples
    /// ```
    /// use rudac::heap::FibonacciHeap;
    ///
    /// let mut fibonacci_heap: FibonacciHeap<usize> = FibonacciHeap::init_min();
    /// fibonacci_heap.push(0);
    ///
    /// fibonacci_heap.clear();
    ///
    /// assert_eq!(fibonacci_heap.size(), 0);
    /// assert_eq!(fibonacci_heap.pop(), None);
    /// ```
    pub fn clear(&mut self) {
        self.children_list.clear();
        self.size = 0;
        self.priority_pointer = None;
    }

    /// Returns true if there are no more items in the heap
    ///
    /// # Examples
    /// ```
    /// use rudac::heap::FibonacciHeap;
    ///
    /// let mut fibonacci_heap: FibonacciHeap<usize> = FibonacciHeap::init_min();
    ///
    /// fibonacci_heap.push(0);
    /// assert_eq!(fibonacci_heap.is_empty(), false);
    ///
    /// fibonacci_heap.pop();
    /// assert_eq!(fibonacci_heap.is_empty(), true);
    /// ```
    pub fn is_empty(&self) -> bool {
        self.size == 0
    }

    /// Returns number of items in heap
    ///
    /// # Examples
    /// ```
    /// use rudac::heap::FibonacciHeap;
    ///
    /// let mut fibonacci_heap: FibonacciHeap<usize> = FibonacciHeap::init_min();
    /// fibonacci_heap.push(0);
    /// fibonacci_heap.push(1);
    ///
    /// assert_eq!(fibonacci_heap.size(), 2);
    /// ```
    pub fn size(&self) -> usize {
        self.size
    }

    /// Returns true if the heap is initialized as a min heap
    ///
    /// # Examples
    /// ```
    /// use rudac::heap::FibonacciHeap;
    ///
    /// let mut fibonacci_heap: FibonacciHeap<usize> = FibonacciHeap::init_min();
    ///
    /// assert_eq!(fibonacci_heap.is_min(), true);
    /// ```
    pub fn is_min(&self) -> bool {
        self.min
    }

    /// Returns true if the heap is initialized as a max heap
    ///
    /// # Examples
    /// ```
    /// use rudac::heap::FibonacciHeap;
    ///
    /// let mut fibonacci_heap: FibonacciHeap<usize> = FibonacciHeap::init_max();
    ///
    /// assert_eq!(fibonacci_heap.is_max(), true);
    /// ```
    pub fn is_max(&self) -> bool {
        !self.is_min()
    }
}
