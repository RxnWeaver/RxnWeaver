package molecule

// InMessage is a message sent to a molecule by an external agent.
//
// An in-message comprises a request, an optional cookie value and a
// dynamic payload.  The meaning of the cookie depends on the request.
// Similarly, the actual contents of the payload depend on the
// request.
//
// Thus, it is highly imperative that other agents that correspond
// with a molecule be aware of what requests molecules understand, and
// what payloads are to be delivered as part of the message.
type InMessage struct {
	Request uint8
	Cookie  uint64
	Payload interface{}
}

// OutMessage is a message sent by a molecule in response to an
// in-message.
//
// An out-message comprises a status (result code), an optional cookie
// value and a dynamic payload.  The meaning of the cookie depends on
// the request.  Similarly, the actual contents of the payload depend
// on the request.
//
// Thus, it is highly imperative that other agents that correspond
// with a molecule be aware of what responses molecules send, and what
// payloads are to be delivered as part of the message.
type OutMessage struct {
	Status  int16
	Cookie  uint64
	Payload interface{}
}
