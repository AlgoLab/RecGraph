use slotmap::{new_key_type, Key, SlotMap};

new_key_type! {
    pub struct NodeKey;
}

#[derive(Copy, Clone)]
struct Node<T> {
    g: T,
    h: T,
    parent: NodeKey,
}

pub struct AlignGraph<T> {
    sm: SlotMap<NodeKey, Node<T>>,
    head: NodeKey,
    tail: NodeKey,
}

impl<T> AlignGraph<T> {
    pub fn new() -> Self {
        Self {
            sm: SlotMap::with_key(),
            head: NodeKey::null(),
            tail: NodeKey::null(),
        }
    }

    pub fn len(&self) -> usize {
        self.sm.len()
    }

    pub fn push_head(&mut self, value: T) -> NodeKey {
        let k = self.sm.insert(Node {
            value,
            prev: NodeKey::null(),
            next: self.head,
        });

        if let Some(old_head) = self.sm.get_mut(self.head) {
            old_head.prev = k;
        } else {
            self.tail = k;
        }
        self.head = k;
        k
    }

    pub fn push_tail(&mut self, value: T) -> NodeKey {
        let k = self.sm.insert(Node {
            value,
            prev: self.tail,
            next: NodeKey::null(),
        });

        if let Some(old_tail) = self.sm.get_mut(self.tail) {
            old_tail.next = k;
        } else {
            self.head = k;
        }
        self.tail = k;
        k
    }

    pub fn pop_head(&mut self) -> Option<T> {
        self.sm.remove(self.head).map(|old_head| {
            self.head = old_head.next;
            old_head.value
        })
    }

    pub fn pop_tail(&mut self) -> Option<T> {
        self.sm.remove(self.tail).map(|old_tail| {
            self.tail = old_tail.prev;
            old_tail.value
        })
    }

    pub fn remove(&mut self, key: NodeKey) -> Option<T> {
        self.sm.remove(key).map(|node| {
            if let Some(prev_node) = self.sm.get_mut(node.prev) {
                prev_node.next = node.next;
            } else {
                self.head = node.next;
            }

            if let Some(next_node) = self.sm.get_mut(node.next) {
                next_node.prev = node.prev;
            } else {
                self.tail = node.prev;
            }

            node.value
        })
    }

    pub fn head(&self) -> NodeKey {
        self.head
    }

    pub fn tail(&self) -> NodeKey {
        self.tail
    }

    pub fn get(&self, key: NodeKey) -> Option<&T> {
        self.sm.get(key).map(|node| &node.value)
    }

    pub fn get_mut(&mut self, key: NodeKey) -> Option<&mut T> {
        self.sm.get_mut(key).map(|node| &mut node.value)
    }
}
