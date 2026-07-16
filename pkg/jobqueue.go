package ggtype

import (
)

type RunQueue struct {
	queue chan struct{}
}

func NewRunQueue(n int) *RunQueue {
	if n < 1 {
		return &RunQueue{queue: make(chan struct{}, 0)}
	}
	return &RunQueue{queue: make(chan struct{}, n)}
}

func (r *RunQueue) Run(f func()) {
	if cap(r.queue) < 1 {
		f()
		return
	}

	r.queue <- struct{}{}
	f()
	<-r.queue
}

func (r *RunQueue) RunErr(f func() error) error {
	var e error
	r.Run(func() {
		e = f()
	})
	return e
}

func (r *RunQueue) Enqueue() {
	if cap(r.queue) < 1 {
		return
	}
	r.queue <- struct{}{}
}

func (r *RunQueue) Dequeue() {
	if cap(r.queue) < 1 {
		return
	}
	<-r.queue
}
